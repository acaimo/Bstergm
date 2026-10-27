#include "ergm_changestat.h"
#include "ergm_storage.h"

#define SA_POS 1u
#define SA_ACT 2u

typedef struct {
  Vertex n;
  Edge npi;
  unsigned char *st;
  Vertex *adeg;
  Vertex **anbr;
  Vertex *anbr_data;
  int *dneg;
  unsigned short *sp[3];
} SignAux;

static inline size_t sa_ut(Vertex a, Vertex b, Vertex n){
  if(a > b){ Vertex t = a; a = b; b = t; }
  size_t i = a - 1, j = b - 1;
  return i * (2 * (size_t) n - i - 1) / 2 + (j - i - 1);
}

static inline unsigned char sa_flag(const SignAux *sa, Vertex a, Vertex b){
  return sa->st[(size_t)(a - 1) * sa->n + (b - 1)];
}

static inline int sa_sign_of(unsigned char f){
  return (f & SA_POS) ? 1 : ((f & SA_ACT) ? -1 : 0);
}

static inline int sa_sign(const SignAux *sa, Vertex a, Vertex b){
  return sa_sign_of(sa_flag(sa, a, b));
}

static inline int sa_match(int p, int u, int v){
  switch(p){
  case 0: return u == 1 && v == 1;
  case 1: return u == -1 && v == -1;
  default: return (u == 1 && v == -1) || (u == -1 && v == 1);
  }
}

static inline int sa_ptype(int u, int v){
  return (u == 1 && v == 1) ? 0 : ((u == -1 && v == -1) ? 1 : 2);
}

#define SA_FOR_SUPPORT(sa, v, k, body)                                    \
  {                                                                       \
    const Vertex *sa_nb_ = (sa)->anbr[(v) - 1];                           \
    const Vertex sa_d_ = (sa)->adeg[(v) - 1];                             \
    for(Vertex sa_q_ = 0; sa_q_ < sa_d_; sa_q_++){                        \
      Vertex k = sa_nb_[sa_q_];                                           \
      body;                                                               \
    }                                                                     \
    if((sa)->npi){                                                        \
      Edge sa_e_;                                                         \
      Vertex k;                                                           \
      STEP_THROUGH_OUTEDGES((v), sa_e_, k){                               \
        if(!(sa_flag((sa), (v), k) & SA_ACT)) body;                       \
      }                                                                   \
      STEP_THROUGH_INEDGES((v), sa_e_, k){                                \
        if(!(sa_flag((sa), (v), k) & SA_ACT)) body;                       \
      }                                                                   \
    }                                                                     \
  }

I_CHANGESTAT_FN(i__signaux){
  if(DIRECTED || BIPARTITE)
    error("Signed network terms require an undirected, unipartite network.");
  Vertex n = N_NODES;
  if(n > 65535) error("Signed network terms support at most 65535 nodes.");

  ALLOC_AUX_STORAGE(1, SignAux, sa);
  sa->n = n;
  sa->npi = 0;
  sa->st = R_Calloc((size_t) n * n, unsigned char);
  sa->adeg = R_Calloc(n, Vertex);
  sa->dneg = R_Calloc(n, int);
  sa->anbr = R_Calloc(n, Vertex *);

  Edge nact = IINPUT_PARAM[0];
  const int *tl = IINPUT_PARAM + 1, *hd = IINPUT_PARAM + 1 + nact;
  unsigned char *keep = R_Calloc(nact ? nact : 1, unsigned char);
  for(Edge e = 0; e < nact; e++){
    Vertex t = tl[e], h = hd[e];
    if(t < 1 || h < 1 || t > n || h > n || t == h) continue;
    unsigned char *f = sa->st + (size_t)(t - 1) * n + (h - 1);
    if(*f & SA_ACT) continue;
    *f |= SA_ACT;
    sa->st[(size_t)(h - 1) * n + (t - 1)] |= SA_ACT;
    sa->adeg[t - 1]++;
    sa->adeg[h - 1]++;
    keep[e] = 1;
  }

  size_t tot = 0;
  for(Vertex v = 0; v < n; v++) tot += sa->adeg[v];
  sa->anbr_data = R_Calloc(tot ? tot : 1, Vertex);
  Vertex *fill = R_Calloc(n, Vertex);
  size_t pos = 0;
  for(Vertex v = 0; v < n; v++){
    sa->anbr[v] = sa->anbr_data + pos;
    pos += sa->adeg[v];
  }
  for(Edge e = 0; e < nact; e++){
    if(!keep[e]) continue;
    Vertex t = tl[e], h = hd[e];
    sa->anbr[t - 1][fill[t - 1]++] = h;
    sa->anbr[h - 1][fill[h - 1]++] = t;
  }
  R_Free(fill);
  R_Free(keep);

  EXEC_THROUGH_NET_EDGES(t, h, e, {
      sa->st[(size_t)(t - 1) * n + (h - 1)] |= SA_POS;
      sa->st[(size_t)(h - 1) * n + (t - 1)] |= SA_POS;
      if(!(sa->st[(size_t)(t - 1) * n + (h - 1)] & SA_ACT)) sa->npi++;
    });

  for(Vertex v = 1; v <= n; v++){
    const Vertex *nb = sa->anbr[v - 1];
    int d = 0;
    for(Vertex q = 0; q < sa->adeg[v - 1]; q++)
      if(!(sa_flag(sa, v, nb[q]) & SA_POS)) d++;
    sa->dneg[v - 1] = d;
  }

  size_t npairs = (size_t) n * (n - 1) / 2;
  for(int p = 0; p < 3; p++) sa->sp[p] = R_Calloc(npairs ? npairs : 1, unsigned short);
  Vertex *buf = R_Calloc(n, Vertex);
  signed char *bsg = R_Calloc(n, signed char);
  for(Vertex m = 1; m <= n; m++){
    Vertex len = 0;
    SA_FOR_SUPPORT(sa, m, k, {
        buf[len] = k;
        bsg[len] = (signed char) sa_sign(sa, m, k);
        len++;
      });
    for(Vertex a = 0; a < len; a++)
      for(Vertex b = a + 1; b < len; b++)
        sa->sp[sa_ptype(bsg[a], bsg[b])][sa_ut(buf[a], buf[b], n)]++;
  }
  R_Free(bsg);
  R_Free(buf);
}

static inline void sa_update_pair(SignAux *sa, Vertex x, Vertex k, int leg, int a, int b){
  size_t idx = sa_ut(x, k, sa->n);
  for(int p = 0; p < 3; p++){
    int d = sa_match(p, b, leg) - sa_match(p, a, leg);
    if(d) sa->sp[p][idx] = (unsigned short)(sa->sp[p][idx] + d);
  }
}

U_CHANGESTAT_FN(u__signaux){
  GET_AUX_STORAGE(SignAux, sa);
  Vertex n = sa->n;
  unsigned char f = sa_flag(sa, tail, head);
  int a = sa_sign_of(f);
  int b = edgestate ? ((f & SA_ACT) ? -1 : 0) : 1;

  if(a != b){
    SA_FOR_SUPPORT(sa, head, k, {
        if(k != tail) sa_update_pair(sa, tail, k, sa_sign(sa, head, k), a, b);
      });
    SA_FOR_SUPPORT(sa, tail, k, {
        if(k != head) sa_update_pair(sa, head, k, sa_sign(sa, tail, k), a, b);
      });
  }

  int dn = (b == -1) - (a == -1);
  sa->dneg[tail - 1] += dn;
  sa->dneg[head - 1] += dn;

  if(!(f & SA_ACT)){
    if(edgestate) sa->npi--;
    else sa->npi++;
  }

  sa->st[(size_t)(tail - 1) * n + (head - 1)] ^= SA_POS;
  sa->st[(size_t)(head - 1) * n + (tail - 1)] ^= SA_POS;
}

F_CHANGESTAT_FN(f__signaux){
  GET_AUX_STORAGE(SignAux, sa);
  if(sa){
    R_Free(sa->st);
    R_Free(sa->adeg);
    R_Free(sa->anbr);
    R_Free(sa->anbr_data);
    R_Free(sa->dneg);
    for(int p = 0; p < 3; p++) R_Free(sa->sp[p]);
    R_Free(sa);
  }
  AUX_STORAGE = NULL;
}

C_CHANGESTAT_FN(c_esp_signed){
  GET_AUX_STORAGE(SignAux, sa);
  const int base = IINPUT_PARAM[0], p = IINPUT_PARAM[1];
  const unsigned int ns = N_CHANGE_STATS;
  const Vertex n = sa->n;
  const size_t L = n > 1 ? n - 1 : 1;
  const double *tab = INPUT_PARAM;
  const unsigned short *sp = sa->sp[p];
  unsigned char f = sa_flag(sa, tail, head);
  int a = sa_sign_of(f);
  int b = edgestate ? ((f & SA_ACT) ? -1 : 0) : 1;

  ZERO_ALL_CHANGESTATS();

  int bs = (b == base) - (a == base);
  if(bs){
    unsigned int s = sp[sa_ut(tail, head, n)];
    for(unsigned int q = 0; q < ns; q++) CHANGE_STAT[q] += bs * tab[q * L + s];
  }

  SA_FOR_SUPPORT(sa, head, k, {
      if(k != tail){
        int leg = sa_sign(sa, head, k);
        int d = sa_match(p, b, leg) - sa_match(p, a, leg);
        if(d && sa_sign(sa, tail, k) == base){
          unsigned int s = sp[sa_ut(tail, k, n)];
          for(unsigned int q = 0; q < ns; q++)
            CHANGE_STAT[q] += tab[q * L + s + d] - tab[q * L + s];
        }
      }
    });
  SA_FOR_SUPPORT(sa, tail, k, {
      if(k != head){
        int leg = sa_sign(sa, tail, k);
        int d = sa_match(p, b, leg) - sa_match(p, a, leg);
        if(d && sa_sign(sa, head, k) == base){
          unsigned int s = sp[sa_ut(head, k, n)];
          for(unsigned int q = 0; q < ns; q++)
            CHANGE_STAT[q] += tab[q * L + s + d] - tab[q * L + s];
        }
      }
    });
}

C_CHANGESTAT_FN(c_degree_neg_signed){
  GET_AUX_STORAGE(SignAux, sa);
  const unsigned int ns = N_CHANGE_STATS;
  const size_t L = sa->n;
  const double *tab = INPUT_PARAM;
  unsigned char f = sa_flag(sa, tail, head);
  int a = sa_sign_of(f);
  int b = edgestate ? ((f & SA_ACT) ? -1 : 0) : 1;
  int dn = (b == -1) - (a == -1);

  ZERO_ALL_CHANGESTATS();
  if(dn){
    int d0 = sa->dneg[tail - 1], d1 = sa->dneg[head - 1];
    for(unsigned int q = 0; q < ns; q++)
      CHANGE_STAT[q] += tab[q * L + d0 + dn] - tab[q * L + d0]
        + tab[q * L + d1 + dn] - tab[q * L + d1];
  }
}

static inline unsigned int sa_tri_count(int npos, int s, unsigned int pp, unsigned int mm, unsigned int pm){
  if(s == 0) return 0;
  int own = s == 1;
  unsigned int r = 0;
  if(npos == own + 2) r += pp;
  if(npos == own + 1) r += pm;
  if(npos == own) r += mm;
  return r;
}

C_CHANGESTAT_FN(c_triangle_signed){
  GET_AUX_STORAGE(SignAux, sa);
  const unsigned int ns = N_CHANGE_STATS;
  size_t idx = sa_ut(tail, head, sa->n);
  unsigned int pp = sa->sp[0][idx], mm = sa->sp[1][idx], pm = sa->sp[2][idx];
  unsigned char f = sa_flag(sa, tail, head);
  int a = sa_sign_of(f);
  int b = edgestate ? ((f & SA_ACT) ? -1 : 0) : 1;

  for(unsigned int q = 0; q < ns; q++){
    int npos = IINPUT_PARAM[q];
    CHANGE_STAT[q] = (double) sa_tri_count(npos, b, pp, mm, pm) - (double) sa_tri_count(npos, a, pp, mm, pm);
  }
}
