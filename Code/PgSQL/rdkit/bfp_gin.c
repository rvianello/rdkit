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

/* gin_bfp_ops opclass options */
typedef struct
{
  int32 vl_len_;		/* varlena header (do not touch directly!) */
  int bands;			  /* number of bands in the minhash signature */
  int rows;			    /* number of rows in each band */
} GinBfpOptions;

#define BANDS_DEFAULT 10
#define BANDS_MIN 8
#define BANDS_MAX 64
#define GET_BANDS()	(PG_HAS_OPCLASS_OPTIONS() ? \
  ((GinBfpOptions *) PG_GET_OPCLASS_OPTIONS())->bands : BANDS_DEFAULT)

#define ROWS_DEFAULT 3
#define ROWS_MIN 1
#define ROWS_MAX 10
#define GET_ROWS()	(PG_HAS_OPCLASS_OPTIONS() ? \
  ((GinBfpOptions *) PG_GET_OPCLASS_OPTIONS())->rows : ROWS_DEFAULT)

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

  for (int32 i = 0; i < nkeys; ++i) {
    if (check[i] == true) {
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

  GinTernaryValue result = GIN_MAYBE;

  for (int32 i = 0; i < nkeys; ++i) {
    if ((check[i] == GIN_TRUE) || (check[i] == GIN_MAYBE)) {
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

  PG_RETURN_VOID();
}
