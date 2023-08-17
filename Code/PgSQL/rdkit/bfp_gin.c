//
//  Copyright (c) 2023, Riccardo Vianello
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of the authors nor the names of their contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <postgres.h>

#include <access/gin.h>
#include <access/reloptions.h>
#include <access/stratnum.h>
#include <fmgr.h>

#include "rdkit.h"
#include "bitstring.h"
#include "guc.h"

/*
 * The implementation is based on indexing the binary fingerprints under locality-sensitive
 * hash keys that are computed from a minhash signature of the original bitstrings.
 *
 * The approach is based on chapter 3 of "The Mining of Massive Datasets" (http://mmds.org/).
 *
 * Given two binary fingerprints, the probability s that their minhash signatures have the
 * same value at any given position, equals their tanimoto similarity. If the signature is
 * divided into b bands or r rows each, the probability p that the signature values match on
 * every row of a given band is s**r. The probability that the two signatures have matching
 * values on each row of k bands then follows a binomial distribution of parameters b (the
 * number of bands) and p = s**r.
 *
 * Searching the index for items that share at least a required number of hash keys with the
 * query is therefore subject an S-shaped probability distribution of selecting the indexed
 * fingerprints based on their increasing tanimoto similarity to the query (the behavior
 * is similar to that of a partial index, including with probability ~1 only the items that are
 * sufficiently similar to the query and excluding those that are too dissimilar).
 *
 * The configured parameters (bands, rows, and min # of matching hash keys) affect the shape
 * of this probability distribution and consequently the selectivity of the index.
 */

/* gin_bfp_ops opclass options */
typedef struct
{
  int32 vl_len_;		     /* varlena header (do not touch directly!) */
  int bands;			       /* number of bands in the minhash signature */
  int rows;			         /* number of rows in each band */
  int required_matching; /* minimum number of matching hash keys (bands) of candidate similar records */
} GinBfpOptions;

#define BANDS_DEFAULT 20
#define BANDS_MIN 8
#define BANDS_MAX 32
#define GET_BANDS()	(PG_HAS_OPCLASS_OPTIONS() ? \
  ((GinBfpOptions *) PG_GET_OPCLASS_OPTIONS())->bands : BANDS_DEFAULT)

#define ROWS_DEFAULT 5
#define ROWS_MIN 1
#define ROWS_MAX 10
#define GET_ROWS()	(PG_HAS_OPCLASS_OPTIONS() ? \
  ((GinBfpOptions *) PG_GET_OPCLASS_OPTIONS())->rows : ROWS_DEFAULT)

#define REQUIRED_MATCHING_DEFAULT 1
#define REQUIRED_MATCHING_MIN 1
#define REQUIRED_MATCHING_MAX 32
#define GET_REQUIRED_MATCHING()	(PG_HAS_OPCLASS_OPTIONS() ? \
  ((GinBfpOptions *) PG_GET_OPCLASS_OPTIONS())->required_matching : REQUIRED_MATCHING_DEFAULT)


static Datum *gin_bfp_extract(Bfp *bfp, int bands, int rows, int32 *nkeys) {

  uint8 *fp = (uint8 *)VARDATA(bfp);
  int32 siglen = BFP_SIGLEN(bfp);

  uint32 *lshkeys = palloc(sizeof(uint32) * bands);
  calcBfpLSHKeys(fp, siglen, bands, rows, lshkeys);
  *nkeys = bands;
  Datum *keys = palloc(sizeof(Datum) * bands);

  for (int i = 0; i < bands; ++i) {
    keys[i] = UInt32GetDatum(lshkeys[i]);
  }

  return keys;
}

PGDLLEXPORT Datum gin_bfp_extract_value(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_extract_value);
Datum gin_bfp_extract_value(PG_FUNCTION_ARGS) {
  Bfp *bfp = PG_GETARG_BFP_P(0);
  int32 *nkeys = (int32 *)PG_GETARG_POINTER(1);

  int bands = GET_BANDS();
  int rows = GET_ROWS();

  PG_RETURN_POINTER(gin_bfp_extract(bfp, bands, rows, nkeys));
}

PGDLLEXPORT Datum gin_bfp_extract_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_extract_query);
Datum gin_bfp_extract_query(PG_FUNCTION_ARGS) {
  Bfp *bfp = PG_GETARG_BFP_P(0);
  int32 *nkeys = (int32 *)PG_GETARG_POINTER(1);
  /* StrategyNumber strategy = PG_GETARG_UINT16(2); */
  /* bool **pmatch = (bool **) PG_GETARG_POINTER(3); */
  /* Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4); */
  /* bool **nullFlags = (bool **) PG_GETARG_POINTER(5); */
  /* int32 *searchMode = (int32 *)PG_GETARG_POINTER(6); */

  int bands = GET_BANDS();
  int rows = GET_ROWS();

  PG_RETURN_POINTER(gin_bfp_extract(bfp, bands, rows, nkeys));
}

PGDLLEXPORT Datum gin_bfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_consistent);
Datum gin_bfp_consistent(PG_FUNCTION_ARGS) {
  bool *check = (bool *)PG_GETARG_POINTER(0);
  /* StrategyNumber strategy = PG_GETARG_UINT16(1); */
  /* Bfp *query = PG_GETARG_BFP_P(2); */
  int32 nkeys = PG_GETARG_INT32(3);
  /* Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4); */
  bool *recheck = (bool *)PG_GETARG_POINTER(5);
  /* Datum * queryKeys = PG_GETARG_POINTER(6); */
  /* bool *nullFlags = (bool *) PG_GETARG_POINTER(7); */

  *recheck = true;

  int required_matching = GET_REQUIRED_MATCHING();
  int match_count = 0;

  for (int32 i = 0; i < nkeys; ++i) {
    if (check[i] == true) {
      ++match_count;
    }
    if (match_count == required_matching) {
      PG_RETURN_BOOL(true);
    }
  }

  PG_RETURN_BOOL(false);
}

PGDLLEXPORT Datum gin_bfp_triconsistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_triconsistent);
Datum gin_bfp_triconsistent(PG_FUNCTION_ARGS) {
  /*

   */
  GinTernaryValue *check = (GinTernaryValue *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = PG_GETARG_UINT16(1);
  /* Bfp *query = PG_GETARG_BFP_P(2); */
  int32 nkeys = PG_GETARG_INT32(3);
  /* Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4); */
  /* Datum * queryKeys = PG_GETARG_POINTER(5); */
  /* bool *nullFlags = (bool *) PG_GETARG_POINTER(6); */

  int required_matching = GET_REQUIRED_MATCHING();
  int match_count = 0;

  for (int32 i = 0; i < nkeys; ++i) {
    if ((check[i] == GIN_TRUE) || (check[i] == GIN_MAYBE)) {
      ++match_count;
    }
    if (match_count == required_matching) {
      PG_RETURN_GIN_TERNARY_VALUE(GIN_MAYBE);
    }
  }

  PG_RETURN_GIN_TERNARY_VALUE(GIN_FALSE);
}


PGDLLEXPORT Datum gin_bfp_compare(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_compare);
Datum gin_bfp_compare(PG_FUNCTION_ARGS)
{
  uint32		a = PG_GETARG_UINT32(0);
  uint32		b = PG_GETARG_UINT32(1);

  if (a > b)
    PG_RETURN_INT32(1);
  else if (a == b)
    PG_RETURN_INT32(0);
  else
    PG_RETURN_INT32(-1);
}


PGDLLEXPORT Datum gin_bfp_options(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_options);
Datum gin_bfp_options(PG_FUNCTION_ARGS)
{
  local_relopts *relopts = (local_relopts *)PG_GETARG_POINTER(0);

  init_local_reloptions(relopts, sizeof(GinBfpOptions));
  add_local_int_reloption(relopts, "bands",
              "number of bands in the minhash signature",
              BANDS_DEFAULT, BANDS_MIN, BANDS_MAX,
              offsetof(GinBfpOptions, bands));
  add_local_int_reloption(relopts, "rows",
              "number of rows in the minhash signature bands",
              ROWS_DEFAULT, ROWS_MIN, ROWS_MAX,
              offsetof(GinBfpOptions, rows));
  add_local_int_reloption(relopts, "required_matching",
              "number of band hash keys that are required to match the query for a record to be considered a candidate similar",
              REQUIRED_MATCHING_DEFAULT, REQUIRED_MATCHING_MIN, REQUIRED_MATCHING_MAX,
              offsetof(GinBfpOptions, required_matching));

  PG_RETURN_VOID();
}
