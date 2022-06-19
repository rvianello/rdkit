//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written
//       permission.
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
#include <fmgr.h>
#include <access/gist.h>
#include <access/skey.h>
#include <utils/memutils.h>

#include <math.h>

#include "rdkit.h"
#include "guc.h"
#include "cache.h"
#include "bitstring.h"

/*
 * Define the compressed Bfp datum representation (GBfp) to be used
 * as entry in the GiST index
 */

typedef struct {
  char vl_len_[4];
  uint16 minWeight;
  uint16 maxWeight;
  uint8 fp[FLEXIBLE_ARRAY_MEMBER]; /* leaf or inner fingerprint data */
} GBfp;

/* compute the memory size of the index entries based on the 
** signature bfp size (used when instantiating uncompressed keys)
*/
#define GBFP_VARSIZE(x) (sizeof(GBfp) + (x))

/* compute the size of the signature bfp size from the size of the index
** entry.
**
** inner entries for keys with trivial union data do not store any
** fingerprints, and their computed siglen is zero.
*/
#define GBFP_SIGLEN(x) (VARSIZE(x) - sizeof(GBfp))

/* check if the arg entry is compressed / has trivial bfp data */
#define GBFP_ALL1(x)	(VARSIZE(x) == sizeof(GBfp))

/* get the pointer of the entry at given pos from a vector of entries
** (used in the union and picksplit methods)
*/
#define GETENTRY(vec,pos) ((GBfp *) DatumGetPointer((vec)->vector[(pos)].key))

/* collect 'key' into 'result' */
static void merge_key(GBfp *result, GBfp *key);

/* estimate the "distance"/"difference" between two keys */
static int keys_distance(GBfp *v1, GBfp *v2);

/* clone 'key' into a new inner key */
static GBfp * copy_key(GBfp *key);

/*
 * Compress method
 *
 * Converts a data item into a format suitable for physical storage in an index page. 
 */

PGDLLEXPORT Datum gbfp_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_compress);
Datum
gbfp_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY *retval;

  /*
  ** On leaf entries, the Bfp is replaced by a GBfp instance, where
  ** the fingerprint is annotated with the precomputed weight (popcount).
  */
  if (entry->leafkey) {
    Bfp *bfp;
    GBfp *gbfp;
    int size, siglen, weight;
  
    bfp = DatumGetBfpP(entry->key);
    
    siglen = BFP_SIGLEN(bfp);
    weight = bitstringWeight(siglen, (uint8 *)VARDATA(bfp));
    
    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
    
    size = GBFP_VARSIZE(siglen);
    
    gbfp = palloc0(size);
    SET_VARSIZE(gbfp, size);
    gbfp->minWeight = gbfp->maxWeight = weight;
    memcpy(gbfp->fp, VARDATA(bfp), siglen);
    
    gistentryinit(*retval, PointerGetDatum(gbfp),
		  entry->rel, entry->page,
		  entry->offset, false);
  }
  /* on inner nodes check if union is trivial and can be "compressed" */
  else {
    GBfp *key, *ckey;
    int size, siglen;
    bool is_all1, gbfp_all1_key;

    key = (GBfp*) DatumGetPointer(entry->key);

    siglen = GBFP_SIGLEN(key);

    gbfp_all1_key = GBFP_ALL1(key);
    is_all1 = gbfp_all1_key ? true : bitstringAllTrue(siglen, key->fp);

    if (is_all1 == gbfp_all1_key) {
      /* nothing to check or do in this case */
      retval = entry;
    }
    else {
      /* the union doesn't need to be explictly stored */
      Assert(is_all1);

      /* create a key of a suitable size. start with an
       * uncompressed estimate, and remove the size for the implicitly
       * stored fp */
      size = GBFP_VARSIZE(siglen) - siglen;

      retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
    
      ckey = palloc0(size);
      SET_VARSIZE(ckey, size);

      ckey->minWeight = key->minWeight;
      ckey->maxWeight = key->maxWeight;

      gistentryinit(*retval, PointerGetDatum(ckey),
		    entry->rel, entry->page,
		    entry->offset, false);
    }
  }
                
  PG_RETURN_POINTER(retval);
}

/*
 * Decompress method
 *
 * Converts the stored representation of a data item into a format that
 * can be manipulated by the other GiST methods in the operator class.
 */

PGDLLEXPORT Datum gbfp_decompress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_decompress);
Datum
gbfp_decompress(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval;
  GBfp *key;

  key = (GBfp *)DatumGetPointer(PG_DETOAST_DATUM(entry->key));

  if (key != (GBfp *)DatumGetPointer(entry->key)) {
    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
    gistentryinit(*retval, PointerGetDatum(key),
		  entry->rel, entry->page,
		  entry->offset, false);
    PG_RETURN_POINTER(retval);
  }

  PG_RETURN_POINTER(entry);
}

/*
 * Union method
 *
 * summarize the information from a set of keys into a single one
 * the result is always an inner key
 */

PGDLLEXPORT Datum gbfp_union(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_union);
Datum
gbfp_union(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  int *size = (int *) PG_GETARG_POINTER(1);
  
  int i;
  GBfp *result, *key;
  
  key = GETENTRY(entryvec, 0);
  result = copy_key(key);
  *size = VARSIZE(result);
  
  for (i = 1; i < entryvec->n; ++i) {
    key = GETENTRY(entryvec, i);
    merge_key(result, key);
  }
  
  PG_RETURN_POINTER(result);
}

/*
 * Same method
 *
 * check if two entries represent the same data
 */

PGDLLEXPORT Datum gbfp_same(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_same);
Datum
gbfp_same(PG_FUNCTION_ARGS)
{
  GBfp *a = (GBfp *)PG_GETARG_POINTER(0);
  GBfp *b = (GBfp *)PG_GETARG_POINTER(1);
  bool *result = (bool *) PG_GETARG_POINTER(2);

  *result =
    (VARSIZE(a) == VARSIZE(b))
    &&
    (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0)
    ;

  PG_RETURN_POINTER(result);
}

/*
 * Penalty method
 *
 * estimate the cost of inserting a new key into an existing one
 * this latter is always an inner key
 */

PGDLLEXPORT Datum gbfp_penalty(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_penalty);
Datum
gbfp_penalty(PG_FUNCTION_ARGS)
{
  GISTENTRY  *origentry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
  float      *penalty = (float *) PG_GETARG_POINTER(2);
  
  Assert(!GIST_LEAF(origentry)); /* always inserting into inner nodes */

  GBfp *origval = (GBfp *) DatumGetPointer(origentry->key);
  int origsiglen = GBFP_SIGLEN(origval);

  GBfp *newval = (GBfp *) DatumGetPointer(newentry->key);  
  int newsiglen = GBFP_SIGLEN(newval);

  if (origsiglen > 0 && newsiglen != origsiglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  /*
   * The computed penalty determines the insertion path for the new
   * entry into the tree.
   *
   * The penalty should be higher
   * - if the weight of the new entry doesn't fall into the weight
   *   interval of the selected subtree
   * - depending on the estimated minimum distance from the entries in
   *   the subtree
   */

  /* weight penalty
   * how far is the weight of the new entry from the subtree's interval?
   */
  float weight_penalty = 0.f;

  if (newval->minWeight < origval->minWeight) {
    weight_penalty += origval->minWeight - newval->minWeight;
  }
  if (newval->maxWeight > origval->maxWeight) {
    weight_penalty += newval->maxWeight - origval->maxWeight;
  }

  /* minimum distance penalty
   * consider the minimum tanimoto distance of the new entry from the fps in this subtree
   */
  float max_similarity = 1.0f;
  if (origsiglen > 0) { // <=> !GBFP_ALL1(origval)
    float nCommon = (float)(bitstringIntersectionWeight(origsiglen, origval->fp, newval->fp));
    if (nCommon > origval->maxWeight) { nCommon = (float)origval->maxWeight; }
    max_similarity = nCommon / newval->minWeight;
  }

  float dist_penalty = 1.f - max_similarity;

  /* overall penalty
   * the distance penalty is in range (0., 1.) and it'd therefore multiplied by the length of
   * the fingerprint so that it's comparable to the weight penalty. when the union fingerprint
   * in the subtree is filled, the similarity criteria becomes less relevant, the distance 
   * penalty goes to zero and the insertion should be basically determined by the weight penalty.
   */
  *penalty = weight_penalty + origsiglen*dist_penalty;

  PG_RETURN_POINTER(penalty);
}

/*
 * Picksplit method
 * 
 * partition a collection of entries into two clusters
 */

typedef struct {
  OffsetNumber pos;
  int32        cost;
} SPLITCOST;

static int
comparecost(const void *va, const void *vb) {
  const SPLITCOST  *a = (const SPLITCOST *) va;
  const SPLITCOST  *b = (const SPLITCOST *) vb;

  if (a->cost == b->cost) {
    return 0;
  }
  
  return (a->cost > b->cost) ? 1 : -1; 
}

PGDLLEXPORT Datum gbfp_picksplit(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_picksplit);
Datum
gbfp_picksplit(PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);

  OffsetNumber j, k;
  OffsetNumber *left, *right;
  OffsetNumber maxoff;
  int32 nbytes;
  
  OffsetNumber seed_1 = 0, seed_2 = 0;
  
  GBfp *datum_l, *datum_r;
  GBfp *gbfpk, *gbfpj;
  int siglen = 0;
  
  int32 size_alpha, size_beta;
  int32 size_waste, waste = -1;
  SPLITCOST *costvector;

  maxoff = entryvec->n - 1;
  nbytes = (maxoff + 2) * sizeof(OffsetNumber);
  
  v->spl_left = (OffsetNumber *) palloc(nbytes);
  left = v->spl_left;
  v->spl_nleft = 0;

  v->spl_right = (OffsetNumber *) palloc(nbytes);
  right = v->spl_right;
  v->spl_nright = 0;

  /* 
  ** select the GBfp pair that are most dissimilar.
  ** 
  ** these will be the seeds of the two subsets.
  */
  for (k = FirstOffsetNumber; k < maxoff; k = OffsetNumberNext(k)) {
    gbfpk = GETENTRY(entryvec, k);
    if (siglen == 0) {
      siglen = GBFP_SIGLEN(gbfpk);
    }
    for (j = OffsetNumberNext(k); j <= maxoff; j = OffsetNumberNext(j)) {
      gbfpj = GETENTRY(entryvec, j);
      size_waste = keys_distance(gbfpk, gbfpj);
      if (size_waste > waste) {
        waste = size_waste;
        seed_1 = k;
        seed_2 = j;
      }
    }
  }

  /*
  ** initialize two empty subsets
  */
  
  if (seed_1 == 0 || seed_2 == 0) {
    /*
    ** all fps were identical and no waste was measured allowing the seeds to
    ** be assigned, so let's just pick the first two
    */
    seed_1 = 1;
    seed_2 = 2;
  }

  /* form initial .. */
  datum_l = copy_key(GETENTRY(entryvec, seed_1));
  datum_r = copy_key(GETENTRY(entryvec, seed_2));

  /* sort before ... */
  costvector = (SPLITCOST *) palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j)) {
    costvector[j - 1].pos = j;
    gbfpj = GETENTRY(entryvec, j);
    size_alpha = keys_distance(datum_l, gbfpj);
    size_beta = keys_distance(datum_r, gbfpj);
    costvector[j - 1].cost = abs(size_alpha - size_beta);
  }
  qsort((void *) costvector, maxoff, sizeof(SPLITCOST), comparecost);

  for (k = 0; k < maxoff; k++) {
    j = costvector[k].pos;
    
    if (j == seed_1) {
      *left++ = j;
      v->spl_nleft++;
      continue;
    }
    else if (j == seed_2) {
      *right++ = j;
      v->spl_nright++;
      continue;
    }
    
    gbfpj = GETENTRY(entryvec, j);
      
    size_alpha = keys_distance(datum_l, gbfpj);
    size_beta = keys_distance(datum_r, gbfpj);
    
    if ((size_alpha < size_beta) ||
        ((size_alpha == size_beta) && (v->spl_nleft < v->spl_nright))) {
      merge_key(datum_l, gbfpj);
      
      *left++ = j;
      v->spl_nleft++;
    }
    else {
      merge_key(datum_r, gbfpj);
      
      *right++ = j;
      v->spl_nright++;
    }
  }
  
  v->spl_ldatum = PointerGetDatum(datum_l);
  v->spl_rdatum = PointerGetDatum(datum_r);
  
  Assert( v->spl_nleft + v->spl_nright == maxoff );

  PG_RETURN_POINTER(v);
}

static bool
gbfp_inner_consistent(BfpSignature *query, GBfp *key, int siglen,
		      StrategyNumber strategy)
{
  bool result;
  double t;
  double nCommon, nDelta;
  double nQuery = (double) query->weight;

  switch(strategy) {
  case RDKitTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    t = getTanimotoLimit();
    /* 
    ** The following inequalities hold
    ** 
    ** Na*t <= Nb <= Na/t
    **
    ** And for the fingerprints in key we have that 
    ** 
    ** minWeight <= Nb <= maxWeight
    **
    ** so if (Na*t > maxWeight) or (Na/t < minWeight) this subtree can be
    ** discarded.
    */
    if ((key->maxWeight < t*nQuery) || (nQuery < t*key->minWeight)) {
      result = false;
    }
    /* The key in the inner node stores the union of the fingerprints 
    ** that populate the child nodes. We use this union to compute an
    ** upper bound to the similarity. If this upper bound is lower than the 
    ** threashold value the subtree may be pruned.
    **
    ** T = Ncommon / (Na + Nb - Ncommon) <= Ncommon / Na
    */
    else {
      nCommon = (double)(GBFP_ALL1(key) ?
        nQuery :
        bitstringIntersectionWeight(siglen, key->fp, query->fp));
      if (nCommon > key->maxWeight) { nCommon = (double)key->maxWeight; }
      result = nCommon >= t*nQuery;
    }
    break;
  case RDKitDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    t = getDiceLimit();
    nCommon = (double)(GBFP_ALL1(key) ?
      nQuery :
      bitstringIntersectionWeight(siglen, key->fp, query->fp));
    if (nCommon > key->maxWeight) { nCommon = (double)key->maxWeight; }
    result = 2.0 * nCommon >= t*(nQuery + nCommon);
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_BOOL(result);
}

static bool
gbfp_leaf_consistent(BfpSignature *query, GBfp *key, int siglen,
		     StrategyNumber strategy)
{
  bool result;
  double t;
  double nCommon;
  double nQuery = (double) query->weight;
  double nKey = (double)key->minWeight;

  switch(strategy) {
  case RDKitTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    t = getTanimotoLimit();
    if ((nKey < t*nQuery) || (nQuery < t*nKey)) {
      result = false;
    }
    else {
      nCommon = (double)bitstringIntersectionWeight(siglen,
						    key->fp, query->fp);
      result = nCommon / (nKey + nQuery - nCommon) >= t;
    }
    break;
  case RDKitDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    t = getDiceLimit();
    nCommon = (double)bitstringIntersectionWeight(siglen,
						  key->fp, query->fp);
    result = 2.0 * nCommon / (nKey + nQuery) >= t;
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_BOOL(result);
}


/*
 * Consistent method
 * 
 * Given an index entry and a query constraint, determines whether the index entry
 * is “consistent” with the query (prunes subtrees that can't match the query, and
 * selects the leaf entries that do)
 */
PGDLLEXPORT Datum gbfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_consistent);
Datum
gbfp_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool *recheck = (bool *) PG_GETARG_POINTER(4);
  bool result;

  GBfp *key;
  BfpSignature *query;
  int siglen, gbfp_siglen;
  
  *recheck = false;
  
  fcinfo->flinfo->fn_extra = searchBfpCache(
					    fcinfo->flinfo->fn_extra,
					    fcinfo->flinfo->fn_mcxt,
					    PG_GETARG_DATUM(1), 
					    NULL, NULL,&query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);
  
  key = (GBfp*) DatumGetPointer(entry->key);
  gbfp_siglen = GBFP_SIGLEN(key);

  if (gbfp_siglen > 0 && siglen != gbfp_siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  result = GIST_LEAF(entry) ?
     gbfp_leaf_consistent(query, key, siglen, strategy) :
     gbfp_inner_consistent(query, key, siglen, strategy)
     ;
  
  PG_RETURN_BOOL(result);
}


static double
gbfp_inner_distance(BfpSignature *query, GBfp *key, int siglen,
		    StrategyNumber strategy)
{
  double nDelta, nQuery, nCommon, similarity;
  
  nQuery = (double)query->weight;
  nCommon = (double)(GBFP_ALL1(key) ?
    nQuery : bitstringIntersectionWeight(siglen, key->fp, query->fp));
  
  if (nCommon > key->maxWeight) {
    nCommon = (double)key->maxWeight;
  }

  switch (strategy) {
  case RDKitOrderByTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    similarity = nCommon / nQuery;
    break;
  case RDKitOrderByDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    similarity =  2.0 * nCommon / (nQuery + nCommon);
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_FLOAT8(1.0 - similarity);
}


static double
gbfp_leaf_distance(BfpSignature *query, GBfp *key, int siglen,
		   StrategyNumber strategy)
{
  double nKey, nQuery, nCommon, similarity;
  
  nKey = (double)key->minWeight;
  nQuery = (double)query->weight;
  nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, query->fp);
    
  switch (strategy) {
  case RDKitOrderByTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    similarity = nCommon / (nKey + nQuery - nCommon);
    break;
  case RDKitOrderByDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    similarity = 2.0 * nCommon / (nKey + nQuery);
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  return 1.0 - similarity;
}


/*
 * Distance method
 * 
 * Given an index entry p and a query value q, this function determines the index
 * entry's “distance” from the query value.
 * For a leaf index entry the result just represents the distance to the index
 * entry; for an internal tree node, the result must be the smallest distance
 * that any child entry could have.
 */
PGDLLEXPORT Datum  gbfp_distance(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_distance);
Datum
gbfp_distance(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  
  GBfp *key = (GBfp *)DatumGetPointer(entry->key);

  BfpSignature *query;
  int siglen, gbfp_siglen;
  double distance;
  
  fcinfo->flinfo->fn_extra = searchBfpCache(fcinfo->flinfo->fn_extra,
					    fcinfo->flinfo->fn_mcxt,
					    PG_GETARG_DATUM(1),
					    NULL, NULL,&query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);
  gbfp_siglen = GBFP_SIGLEN(key);

  if (gbfp_siglen > 0 && siglen != gbfp_siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  distance = GIST_LEAF(entry) ?
    gbfp_leaf_distance(query, key, siglen, strategy) :
    gbfp_inner_distance(query, key, siglen, strategy)
    ;
  
  PG_RETURN_FLOAT8(distance);
}

/*
 * Fetch method
 * 
 * Converts the compressed index representation of a data item into the
 * original data type, for index-only scans.
 */
PGDLLEXPORT Datum  gbfp_fetch(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_fetch);
Datum
gbfp_fetch(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GBfp *gbfp = (GBfp *) DatumGetPointer(PG_DETOAST_DATUM(entry->key));

  int siglen, size;
  Bfp *bfp;
  GISTENTRY  *retval;

  Assert(GIST_LEAF(entry));

  siglen = GBFP_SIGLEN(gbfp);
  
  size = VARHDRSZ + siglen;
    
  bfp = palloc(size);
  SET_VARSIZE(bfp, size);
  memcpy(VARDATA(bfp), gbfp->fp, siglen);
    
  retval = palloc(sizeof(GISTENTRY));
  
  gistentryinit(*retval, PointerGetDatum(bfp),
		entry->rel, entry->page, entry->offset, false);
  
  PG_RETURN_POINTER(retval);
}



/* utility functions */

/*
 * merge_key
 * used in the union and picksplit methods
 * note: the result key must always be an inner key
 */
static void
merge_key(GBfp *result, GBfp *key)
{
  int i, siglen, key_siglen;
  uint8 *fp, *fp_end;

  siglen = GBFP_SIGLEN(result);

  /*
   * we want to update the summary information in the result
   * (inner key) with the data from the new key.
   */
  key_siglen = GBFP_SIGLEN(key);
  if (siglen > 0 && key_siglen > 0 && key_siglen != siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  /* update the weight interval */
  if (key->minWeight < result->minWeight) {
    result->minWeight = key->minWeight;
  }
  if (key->maxWeight > result->maxWeight) {
    result->maxWeight = key->maxWeight;
  }

  /*
   * merging the union fingerprint data must consider that both
   * the result and new key may have trivial data
   */
  if (GBFP_ALL1(result)) {
    /* no need to merge the union fp from key */
  }
  else if (GBFP_ALL1(key)) {
    /* fill the fp union data in result with 1s */
    fp = result->fp;
    fp_end = fp + siglen;
    while (fp < fp_end) {
      *fp++ = 0xff;
    }
  }
  else {
    /* merge the fp union data from key into result */
    bitstringUnion(siglen, result->fp, key->fp);
  }
}

/*
 * keys_distance
 * used in the picksplit methods
 */
static int
keys_distance(GBfp *v1, GBfp *v2)
{
  int32 minw1, maxw1, minw2, maxw2;
  int distance;
  bool is_all1_v1, is_all1_v2;
  
  int siglen = GBFP_SIGLEN(v1);
  int v2_siglen = GBFP_SIGLEN(v2);
  if (siglen > 0 && v2_siglen > 0 && v2_siglen != siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  /*
   * The computed distance includes contributions coming from
   * - the distance between the weight intervals
   * - the distance between the union bfps
   */

  /* weight distance */
  int dwmin = abs(v1->minWeight - v2->minWeight);
  int dwmax = abs(v1->maxWeight - v2->maxWeight);

  distance = dwmin;
  if (dwmax > distance) {
    distance = dwmax;
  }
  
#if 0
  minw1 = v1->minWeight;
  maxw1 = v1->maxWeight;
    
  minw2 = v2->minWeight;
  maxw2 = v2->maxWeight;

  distance = abs(minw1 - minw2) + abs(maxw1 - maxw2);

  is_all1_v1 = GBFP_ALL1(v1);
  is_all1_v2 = GBFP_ALL1(v2);

  /* union bfp distance */

  if (is_all1_v1 && is_all1_v2) {
    /* both keys have trivial union fp data */
    /* distance += 0; */
  }
  else if (is_all1_v1) {
    /* v1 has trivial union fp data */
    distance += 8*siglen - bitstringWeight(siglen, v2->fp);
  }
  else if (is_all1_v2) {
    /* v2 has trivial union fp data */
    distance += 8*siglen - bitstringWeight(siglen, v1->fp);
  }
  else {
    /* union bfp data available for both keys */
    distance += bitstringHemDistance(siglen, v1->fp, v2->fp);
  }
#endif
  return distance;
}


/*
 * copy_key
 * used in the union and picksplit methods
 * Note: the result of the copy is always an inner key
 */
static GBfp *
copy_key(GBfp *key)
{
  int size = VARSIZE(key);
  GBfp * result = palloc(size);
  memcpy(result, key, size);
  return result;
}

