#include <stdint.h>
#include <openssl/ec.h>
#include <openssl/bn.h>
#include <openssl/obj_mac.h>
#include <openssl/sha.h>
#include <time.h>
#include <pari/pari.h>

#include "param.h"


/********************************/ 
/*								*/
/* 		 PARI/GP CODE FOR  		*/
/*			3SQR DECOMP			*/
/********************************/

/*
GP;install("init_threesquare","v","init_threesquare","/Users/dahmungoudarzi/git_article/range_proofs/threesquare.gp.so");
GP;install("threesquareRSd","D0,G,p","threesquareRSd","/Users/dahmungoudarzi/git_article/range_proofs/threesquare.gp.so");
*/
void init_threesquare(void);
GEN threesquareRSd(GEN n, long prec);
/*End of prototype*/

void
init_threesquare(void)	  /* void */
{
	return;
}

GEN
threesquareRSd(GEN n, long prec)
{
  GEN nn = gen_0, fourpower = gen_0, introot = gen_0, success = gen_0, ctr = gen_0, u = gen_0, x = gen_0, irt = gen_0, QF = gen_0, p = pol_x(fetch_user_var("p"));
  nn = gcopy(n);
  fourpower = gen_0;
  success = gen_0;
  ctr = gen_0;
  QF = Qfb0(gen_1, gen_0, gen_1, NULL, prec);
  while (gequal0(gmodgs(nn, 4)))
  {
    nn = gdivgs(nn, 4);
    fourpower = gaddgs(fourpower, 1);
  }
  if (gequalgs(gmodgs(nn, 8), 7))
    pari_err(e_MISC, "not a sum of three squares!");
  irt = sqrtint(nn);
  if (gequal(gsqr(irt), nn))
  {
    GEN p1;	  /* vec */
    p1 = cgetg(4, t_VEC);
    gel(p1, 1) = gcopy(irt);
    gel(p1, 2) = gen_0;
    gel(p1, 3) = gen_0;
    return gmul(gpow(gen_2, fourpower, prec), p1);
  }
  if (gequalgs(gmodgs(nn, 4), 2))
  {
    introot = gfloor(gdivgs(subis(sqrtint(nn), 1), 2));
    while (gequal0(success))
    {
      x = gaddgs(gmulsg(2, gsub(introot, ctr)), 1);
      p = gsub(nn, gsqr(x));
      if ((gcmpgs(p, 10000) < 0) || !gequal0(gispseudoprime(p, 0)))
      {
        u = qfbsolve(QF, p, 0);
        if (gequal(gnorml2(u), p))
        {
          GEN p2;	  /* vec */
          p2 = cgetg(4, t_VEC);
          gel(p2, 1) = gcopy(x);
          gel(p2, 2) = gcopy(gel(u, 1));
          gel(p2, 3) = gcopy(gel(u, 2));
          return gmul(gpow(gen_2, fourpower, prec), p2);
        }
      }
      ctr = gaddgs(ctr, 1);
    }
  }
  if (gequal1(gmodgs(nn, 4)))
  {
    introot = gfloor(gdivgs(sqrtint(nn), 2));
    while (gequal0(success))
    {
      x = gmulsg(2, gsub(introot, ctr));
      p = gsub(nn, gsqr(x));
      if (!gequal0(gispseudoprime(p, 0)))
      {
        GEN p3;	  /* vec */
        u = qfbsolve(QF, p, 0);
        p3 = cgetg(4, t_VEC);
        gel(p3, 1) = gcopy(x);
        gel(p3, 2) = gcopy(gel(u, 1));
        gel(p3, 3) = gcopy(gel(u, 2));
        return gmul(gpow(gen_2, fourpower, prec), p3);
      }
      ctr = gaddgs(ctr, 1);
    }
  }
  if (gequalgs(gmodgs(nn, 8), 3))
  {
    introot = gfloor(gdivgs(subis(sqrtint(nn), 1), 2));
    while (gequal0(success))
    {
      x = gaddgs(gmulsg(2, gsub(introot, ctr)), 1);
      p = gdivgs(gsub(nn, gsqr(x)), 2);
      if (!gequal0(gispseudoprime(p, 0)))
      {
        GEN p4;	  /* vec */
        u = qfbsolve(QF, p, 0);
        p4 = cgetg(4, t_VEC);
        gel(p4, 1) = gcopy(x);
        gel(p4, 2) = gadd(gel(u, 1), gel(u, 2));
        gel(p4, 3) = gsub(gel(u, 1), gel(u, 2));
        return gmul(gpow(gen_2, fourpower, prec), p4);
      }
      ctr = gaddgs(ctr, 1);
    }
  }
  return gen_0;
}


/********************************/ 
/*								*/
/* 		 Batch-PoSO IMPL  		*/
/*								*/
/********************************/

int protocol_poso(
	EC_GROUP *group_Gcom,
	EC_GROUP *group_Hcom,
    EC_POINT *commitment_x,
    BIGNUM *batch_x[N],
	BIGNUM *r_x,
	BIGNUM *scalar_gamma,
	BIGNUM *batch_gamma[4][R][N],
	EC_POINT *G_i[N], 
	EC_POINT *G_ij[N][3],
	EC_POINT *Gtilde[R],
	EC_POINT *H_i[N],
	double timing[2]
) {

    BN_CTX *bn_ctx = BN_CTX_new();

  	BIGNUM *range_B = BN_new();
	BN_set_word(range_B, 0xFFFFFFFF);

	BIGNUM *pmod = BN_new(); 
	EC_GROUP_get_order(group_Gcom, pmod, NULL);


	/* PHASE 1 */

	struct timeval start, end;
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	start = usage.ru_utime;

	
	BIGNUM *batch_y[N][3];

	BIGNUM *sqr_in = BN_new();
	BIGNUM *bn_b_min_x = BN_new();
	BIGNUM *bn_4_x = BN_new();
	BIGNUM *bn_prod = BN_new();
	BIGNUM *int_4 = BN_new(), *int_1 = BN_new(), *int_2 = BN_new(), *int_8 = BN_new();
	BN_set_word(int_1, 1); BN_set_word(int_2, 2); BN_set_word(int_4, 4); BN_set_word(int_8, 8);

	// allocates HEAP to pari.
  	pari_init(1000000000000, 2);
	GEN sum_y_2;
    GEN a;
	
	for (size_t i=0; i<N; i++) {
		for (size_t j=0; j<3; j++) {
			batch_y[i][j] = BN_new();
		}
		BN_sub(bn_b_min_x, range_B, batch_x[i]);			// B - xi
		BN_mul(bn_4_x, batch_x[i], int_4, bn_ctx);			// 4xi
		BN_mul(bn_prod, bn_4_x, bn_b_min_x, bn_ctx);		// 4xi (B-xi) 
		BN_add(sqr_in, bn_prod, int_1);						// 4xi (B-xi) + 1
		char *sqrin_str = BN_bn2dec(sqr_in); 
		sum_y_2 = readseq(sqrin_str);				// y^2 = 4xi (B-xi) +1
		a = threesquareRSd(sum_y_2, B);					// a0,a1,a2 <- 3 sqr
		char *y0 = GENtostr(gel(a,1));
		char *y1 = GENtostr(gel(a,2));
		char *y2 = GENtostr(gel(a,3));
		BN_dec2bn(&batch_y[i][0], y0);		
		BN_dec2bn(&batch_y[i][1], y1);
		BN_dec2bn(&batch_y[i][2], y2);
		BN_set_negative(batch_y[i][0], 0);
		BN_set_negative(batch_y[i][1], 0);
		BN_set_negative(batch_y[i][2], 0);

		OPENSSL_free(sqrin_str);
		free(y0);
		free(y1);
		free(y2);
	}
	pari_close();
		
	//// Set commitment for the y's.
	// draw the mu
	BIGNUM *mu[R];
	for (size_t k=0; k<R; k++) {
		mu[k] = BN_new();
		BN_rand(mu[k], 4 + N + B + Gamma + L, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);  // mu[k] <- [0, 4NBGammaL]
	}
    // draw random_y.
	BIGNUM *r_y = BN_new();
	BN_rand(r_y, S, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);					// r_y <- [0, S]

	// compute commitments for the y's
	EC_POINT *commitment_y = EC_POINT_new(group_Gcom);
	EC_POINTs_mul(group_Gcom, commitment_y, r_y, N*3, (const EC_POINT **) G_ij, (const BIGNUM **) batch_y, bn_ctx);     // C_y = r_yG_0 + sum_i^Nsum_j^3 y_i,jG_i,j
	// Add the sum of mu_k Gtilde_k
	EC_POINT *point_tmp_gtilde = EC_POINT_new(group_Gcom);
	for (size_t k=0; k<R; k++) {
		EC_POINT_mul(group_Gcom, point_tmp_gtilde, NULL, Gtilde[k], mu[k], bn_ctx);
		EC_POINT_add(group_Gcom, commitment_y, commitment_y, point_tmp_gtilde, bn_ctx);
	}

	//// Masking and range check
	// Masking
	BIGNUM *zeta[R];
	BIGNUM *zeta_tmp = BN_new();
	int range_zeta;
	for (size_t k=0; k<R; k++) {
		zeta[k] = BN_dup(mu[k]);
		for (size_t i=0; i<N; i++) {
			BN_mul(zeta_tmp, batch_x[i], batch_gamma[3][k][i], bn_ctx);
			BN_add(zeta[k], zeta[k], zeta_tmp);
			for (size_t j=0; j<3; j++) {
				BN_mul(zeta_tmp, batch_y[i][j], batch_gamma[j][k][i], bn_ctx);
				BN_add(zeta[k], zeta[k], zeta_tmp);
			}
		}
		// Range Check
		range_zeta = BN_num_bits(zeta[k]);
		if ((range_zeta < (4 + N + B + Gamma)) && ((4 + N + B + Gamma + L) > range_zeta)) {
			// Abort
			return 1;
		}
	}
	
	/* PHASE 2 */

	//// Commitment to random masks
    // draw random from [0, p-1]
	BIGNUM *rtilde_x = BN_new(), *rtilde_y = BN_new();
	BN_rand_range(rtilde_x, pmod);
	BN_rand_range(rtilde_y, pmod);

	// draw random from [0, BOmegaL]
	BIGNUM *xtilde_i[N], *ytilde_ij[N][3];
	for (size_t i=0; i<N; i++) {
		xtilde_i[i]= BN_new();
		BN_rand(xtilde_i[i], B + Omega + L, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
		for (size_t j=0; j<3; j++) {
			ytilde_ij[i][j]= BN_new();
			BN_rand(ytilde_ij[i][j], B + Omega + L, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
		}
	}

	BIGNUM *mutilde[R];
	for (size_t k=0; k<R; k++) {
		mutilde[k] = BN_new();
		BN_rand_range(mutilde[k], pmod);
	}
	
	//// Compute random masks D_x, D_y, and d_k
	// Dx = rtilde_xG_0 + sum_i^N x_iG_i
    EC_POINT *mask_Dx = EC_POINT_new(group_Gcom);
	EC_POINTs_mul(group_Gcom, mask_Dx, rtilde_x, N, (const EC_POINT **) G_i, (const BIGNUM **) xtilde_i, bn_ctx);


	// Dy = rtilde_yG_0 + sum_i^N sum_j^3 ytilde_ijG_ij + sum_k^R mutilde_kGtilde_k
	EC_POINT *mask_Dy = EC_POINT_new(group_Gcom);
	EC_POINTs_mul(group_Gcom, mask_Dy, rtilde_y, N*3, (const EC_POINT **) (const BIGNUM **) G_ij, (const BIGNUM **) ytilde_ij, bn_ctx);
	// Add the sum of mutilde_k Gtilde_k
	for (size_t k=0; k<R; k++) {
		EC_POINT_mul(group_Gcom, point_tmp_gtilde, NULL, Gtilde[k], mutilde[k], bn_ctx);
		EC_POINT_add(group_Gcom, mask_Dy, mask_Dy, point_tmp_gtilde, bn_ctx);
	}

	// d_k = sum_i^Nsum_j^3 ytilde_ijgamma_ijk + mutilde_k
	BIGNUM *dk[R];
	BIGNUM *dk_tmp = BN_new();
	for (size_t k=0; k<R; k++) {
		dk[k] = BN_dup(mutilde[k]);
		for (size_t i=0; i<N; i++) {
			BN_mul(dk_tmp, xtilde_i[i], batch_gamma[3][k][i], bn_ctx);
			BN_add(dk[k], dk[k], dk_tmp);
			for (size_t j=0; j<3; j++) {
				BN_mul(dk_tmp, ytilde_ij[i][j], batch_gamma[j][k][i], bn_ctx);
				BN_add(dk[k], dk[k], dk_tmp);
			}
		}
		BN_mod(dk[k], dk[k], pmod, bn_ctx);
	}

	//// Commitment to decomposition coefficients
	BIGNUM *alpha_one[N];
    BIGNUM *alpha_zero[N];
    BIGNUM *r_star = BN_new(), *rtilde_star = BN_new(); 

	BIGNUM *mB = BN_new();
	BIGNUM *x8 = BN_new();
	BIGNUM *xm = BN_new();
	BIGNUM *ym = BN_new();
	BIGNUM *ym0 = BN_new();
	BIGNUM *ym1 = BN_new();
	BIGNUM *ym2 = BN_new();

	BIGNUM *mki2 = BN_new();
	BIGNUM *mkij2_0 = BN_new();
	BIGNUM *mkij2_1 = BN_new();
	BIGNUM *mkij2_2 = BN_new();

	BN_rand(r_star, S, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
	BN_rand_range(rtilde_star, pmod);
	for (size_t i=0; i<N; i++) {
		// alpha_one computation
		BN_mul(mB, xtilde_i[i], range_B, bn_ctx);
		BN_mul(mB, mB, int_4, bn_ctx);
		BN_mul(x8, batch_x[i], int_8, bn_ctx);
		BN_mul(xm, x8, xtilde_i[i], bn_ctx);
		BN_mul(ym0, batch_y[i][0], ytilde_ij[i][0], bn_ctx);
		BN_mul(ym1, batch_y[i][1], ytilde_ij[i][1], bn_ctx);
		BN_mul(ym2, batch_y[i][2], ytilde_ij[i][2], bn_ctx);
		BN_add(ym, ym0, ym1);
		BN_add(ym, ym, ym2);
		BN_mul(ym, ym, int_2, bn_ctx);

		alpha_one[i] = BN_new();
		BN_sub(alpha_one[i], mB, xm);
		BN_sub(alpha_one[i], alpha_one[i], ym);

		// alpha_zero computation
		BN_sqr(mki2, xtilde_i[i], bn_ctx);
		BN_mul(mki2, mki2, int_4, bn_ctx);
		BN_sqr(mkij2_0, ytilde_ij[i][0], bn_ctx);
		BN_sqr(mkij2_1, ytilde_ij[i][1], bn_ctx);
		BN_sqr(mkij2_2, ytilde_ij[i][2], bn_ctx);

		alpha_zero[i] = BN_new();
		BN_add(alpha_zero[i], mkij2_0, mkij2_1);
		BN_add(alpha_zero[i], alpha_zero[i], mkij2_2);
		BN_add(alpha_zero[i], alpha_zero[i], mki2);
		BN_set_negative(alpha_zero[i], 1);

	}

		
	// Cstar and Dstar computation
	EC_POINT *C_star= EC_POINT_new(group_Hcom), *D_star = EC_POINT_new(group_Hcom);

	EC_POINTs_mul(group_Hcom, C_star, r_star, N, (const EC_POINT **) H_i, (const BIGNUM **) alpha_one, bn_ctx);
	EC_POINTs_mul(group_Hcom, D_star, rtilde_star, N, (const EC_POINT **) H_i, (const BIGNUM **) alpha_zero, bn_ctx);


	// Compute Delta = Hash(Dx, Dy, D*) 
	char *char_dk[R];
	size_t size_dx, size_dy, size_dstar, size_dk = 0;
	size_t size_concatD = 0;

	
	char *char_Dx = EC_POINT_point2hex(group_Gcom, mask_Dx, POINT_CONVERSION_COMPRESSED, bn_ctx);
	char *char_Dy = EC_POINT_point2hex(group_Gcom, mask_Dy, POINT_CONVERSION_COMPRESSED, bn_ctx);
	char *char_Dstar = EC_POINT_point2hex(group_Hcom, D_star, POINT_CONVERSION_COMPRESSED, bn_ctx);
	for (size_t k=0; k<R; k++) {
		char_dk[k] = BN_bn2hex(dk[k]);
		size_dk += strlen(char_dk[k]);
	}
	size_dx = strlen(char_Dx);
	size_dy = strlen(char_Dy);
	size_dstar = strlen(char_Dstar);
	size_concatD = size_dx + size_dy + size_dstar + size_dk;

	char concat_D[size_concatD];
	strcpy(concat_D, char_Dx);
	strcat(concat_D, char_Dy);
	strcat(concat_D, char_Dstar);
	for (size_t k=0; k<R; k++) {
		strcat(concat_D, char_dk[k]);
	}
	unsigned char delta[32];
	SHA256((const unsigned char*) concat_D, size_concatD, delta);

	//// Response
	BIGNUM *gx = BN_new();
	BIGNUM *zi[N], *zij[N][3];
	BIGNUM *tx, *ty, *tstar, *tau[R];

	int range_zi, range_zij;

	for (size_t i=0; i<N; i++) {
		// mask_zi
		zi[i] = BN_new();
		BN_mul(gx, scalar_gamma, batch_x[i], bn_ctx);
		BN_add(zi[i], gx, xtilde_i[i]);
		BN_clear(gx);
		range_zi = BN_num_bits(zi[i]);
		if ((range_zi < (B + Omega)) && ((B + Omega + L) > range_zi)) {
			// Abort
			return 1;
		}
		for (size_t j=0; j<3; j++) {
			// mask_zij
			zij[i][j] = BN_new();
			BN_mul(gx, scalar_gamma, batch_y[i][j], bn_ctx);
			BN_add(zij[i][j], gx, ytilde_ij[i][j]);
			range_zij = BN_num_bits(zij[i][j]);
			BN_clear(gx);
			if ((range_zij < (B + Omega)) && ((B + Omega + L) > range_zij)) {
				// Abort
				return 1;
			}
		}
	}

	// mask_tx
	tx = BN_new();
	BN_mul(gx, scalar_gamma, r_x, bn_ctx);
	BN_add(tx, gx, rtilde_x);
	BN_clear(gx);
	BN_mod(tx, tx, pmod, bn_ctx);

	// mask_ty
	ty = BN_new();
	BN_mul(gx, scalar_gamma, r_y, bn_ctx);
	BN_add(ty, gx, rtilde_y);
	BN_clear(gx);
	BN_mod(ty, ty, pmod, bn_ctx);

	//mask_tstar
	tstar = BN_new();
	BN_mul(gx, scalar_gamma, r_star, bn_ctx);
	BN_add(tstar, gx, rtilde_star);
	BN_clear(gx);
	BN_mod(tstar, tstar, pmod, bn_ctx);

	//mask tau
	for (size_t k=0; k<R; k++) {
		tau[k] = BN_new();
		BN_mul(gx, scalar_gamma, mu[k], bn_ctx);
		BN_add(tau[k], gx, mutilde[k]);
		BN_mod(tau[k], tau[k], pmod, bn_ctx);
		BN_clear(gx);
	}


	getrusage(RUSAGE_SELF, &usage);
	end = usage.ru_utime;
	double timing_prove = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e6;
	
	//// Verify
	start = end;

	EC_POINT *Fx, *Fy, *Fstar;
	EC_POINT *tmpG = EC_POINT_new(group_Gcom);
	EC_POINT *tmpH = EC_POINT_new(group_Hcom);

	BIGNUM *f_star[N], *fk[R];
	BIGNUM *m_gamma = BN_new();

	BIGNUM *gamma_sqr = BN_new();
	BIGNUM *zij_sqr0 = BN_new();
	BIGNUM *zij_sqr1 = BN_new();
	BIGNUM *zij_sqr2 = BN_new();
	BIGNUM *gammaB = BN_new();
	BIGNUM *zi_4 = BN_new();


	BN_copy(m_gamma, scalar_gamma);
	BN_set_negative(m_gamma, 1);
	Fx = EC_POINT_new(group_Gcom);
	Fy = EC_POINT_new(group_Gcom);

	// Fx = -gammaCx + txG0 + sum_i^N ziGi
	EC_POINTs_mul(group_Gcom, tmpG, tx, N, (const EC_POINT **) G_i, (const BIGNUM **) zi, bn_ctx);		
	EC_POINT_mul(group_Gcom, Fx, NULL, commitment_x, m_gamma, bn_ctx);
	EC_POINT_add(group_Gcom, Fx, Fx, tmpG, bn_ctx);

	// Fy = -gammaCy + tyG0 + sum_i^Nsum_j^3 zijGij + sum_k^RtaukGtilte_k
	EC_POINTs_mul(group_Gcom, tmpG, ty, N*3, (const EC_POINT **) G_ij, (const BIGNUM **) zij, bn_ctx);
	EC_POINT_mul(group_Gcom, Fy, NULL, commitment_y, m_gamma, bn_ctx);
	EC_POINT_add(group_Gcom, Fy, Fy, tmpG, bn_ctx);
	for (size_t k=0; k<R; k++) {
		EC_POINT_mul(group_Gcom, tmpG, NULL, Gtilde[k], tau[k], bn_ctx);
		EC_POINT_add(group_Gcom, Fy, Fy, tmpG, bn_ctx);
	}

	BIGNUM *fk_tmp = BN_new();
	for (size_t k=0; k<R; k++) {
		fk[k] = BN_new();
		BN_mul(fk_tmp, m_gamma, zeta[k], bn_ctx);
		BN_add(fk[k], fk_tmp, tau[k]);
		for (size_t i=0; i<N; i++) {
			BN_mul(fk_tmp, zi[i], batch_gamma[3][k][i], bn_ctx);
			BN_add(fk[k], fk[k], fk_tmp);
			for (size_t j=0; j<3; j++) {
				BN_mul(fk_tmp, zij[i][j], batch_gamma[j][k][i], bn_ctx);
				BN_add(fk[k], fk[k], fk_tmp);
			}
		}
		BN_mod(fk[k], fk[k], pmod, bn_ctx);
	}

	for (size_t i=0; i<N; i++) {
		f_star[i] = BN_new();
		BN_sqr(gamma_sqr, scalar_gamma, bn_ctx);
		BN_sqr(zij_sqr0, zij[i][0], bn_ctx);
		BN_sqr(zij_sqr1, zij[i][1], bn_ctx);
		BN_sqr(zij_sqr2, zij[i][2], bn_ctx);
		BN_mul(zi_4, zi[i], int_4, bn_ctx);
		BN_mul(gammaB, scalar_gamma, range_B, bn_ctx);
		BN_sub(gammaB, gammaB, zi[i]);
		BN_mul(f_star[i], zi_4, gammaB, bn_ctx);
		BN_add(f_star[i], f_star[i], gamma_sqr);
		BN_sub(f_star[i], f_star[i], zij_sqr0);
		BN_sub(f_star[i], f_star[i], zij_sqr1);
		BN_sub(f_star[i], f_star[i], zij_sqr2);
		BN_mod(f_star[i], f_star[i], pmod, bn_ctx);
	}

	Fstar = EC_POINT_new(group_Hcom);
	EC_POINTs_mul(group_Hcom, tmpH, tstar, N, (const EC_POINT **) H_i, (const BIGNUM **) f_star, bn_ctx);
	EC_POINT_mul(group_Hcom, Fstar, NULL, C_star, m_gamma, bn_ctx);
	EC_POINT_add(group_Hcom, Fstar, Fstar, tmpH, bn_ctx);

	//// Verification checks

	// Compute Delta = Hash(Fx, Fy, F*) 
	char *char_fk[R];
	size_t size_fx, size_fy, size_fstar, size_fk=0;
	size_t size_concatF = 0;

	char * char_Fx = EC_POINT_point2hex(group_Gcom, Fx, POINT_CONVERSION_COMPRESSED, bn_ctx);
	char * char_Fy = EC_POINT_point2hex(group_Gcom, Fy, POINT_CONVERSION_COMPRESSED, bn_ctx);
	char * char_Fstar = EC_POINT_point2hex(group_Hcom, Fstar, POINT_CONVERSION_COMPRESSED, bn_ctx);
	for (size_t k=0; k<R; k++) {
		char_fk[k] = BN_bn2hex(fk[k]);
		size_fk += strlen(char_fk[k]);
	}
	size_fx = strlen(char_Fx);
	size_fy = strlen(char_Fy);
	size_fstar = strlen(char_Fstar);
	size_concatF += size_fx + size_fy + size_fstar + size_fk;

	char concat_F[size_concatF];
	strcpy(concat_F, char_Fx);
	strcat(concat_F, char_Fy);
	strcat(concat_F, char_Fstar);
	for (size_t k=0; k<R; k++) {
		strcat(concat_F, char_fk[k]);
	}
	unsigned char delta_verif[32];
	SHA256((const unsigned char*) concat_F, size_concatF, delta_verif);
	
	int verification = memcmp(delta, delta_verif, 32);

	if (verification != 0) {
		return 3;
	}

	getrusage(RUSAGE_SELF, &usage);
	end = usage.ru_utime;
	double timing_verif = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e6;


	timing[0] += timing_prove;
	timing[1] += timing_verif;

	// Clean up
	BN_clear_free(range_B);
	BN_clear_free(pmod);
	BN_clear_free(sqr_in);
	BN_clear_free(bn_b_min_x);
	BN_clear_free(bn_4_x);
	BN_clear_free(bn_prod);
	BN_clear_free(int_8);
	BN_clear_free(int_4);
	BN_clear_free(int_2);
	BN_clear_free(int_1);
	BN_clear_free(r_y);
	BN_clear_free(zeta_tmp);
	BN_clear_free(rtilde_x);
	BN_clear_free(rtilde_y);
	BN_clear_free(dk_tmp);
	BN_clear_free(r_star);
	BN_clear_free(rtilde_star);
	BN_clear_free(mB);
	BN_clear_free(x8);
	BN_clear_free(xm);
	BN_clear_free(ym);
	BN_clear_free(ym0);
	BN_clear_free(ym1);
	BN_clear_free(ym2);
	BN_clear_free(mki2);
	BN_clear_free(mkij2_0);
	BN_clear_free(mkij2_1);
	BN_clear_free(mkij2_2);
	BN_clear_free(gx);
	BN_clear_free(tx);
	BN_clear_free(ty);
	BN_clear_free(tstar);
	BN_clear_free(m_gamma);
	BN_clear_free(gamma_sqr);
	BN_clear_free(zij_sqr0);
	BN_clear_free(zij_sqr1);
	BN_clear_free(zij_sqr2);
	BN_clear_free(gammaB);
	BN_clear_free(zi_4);
	BN_clear_free(fk_tmp);

	EC_POINT_free(commitment_y);
	EC_POINT_free(point_tmp_gtilde);
	EC_POINT_free(mask_Dx);
	EC_POINT_free(mask_Dy);
	EC_POINT_free(C_star);
	EC_POINT_free(D_star);
	EC_POINT_free(tmpG);
	EC_POINT_free(tmpH);
	EC_POINT_free(Fx);
	EC_POINT_free(Fy);
	EC_POINT_free(Fstar);

	for (size_t i=0; i<N; i++) {
		BN_clear_free(xtilde_i[i]);
		BN_clear_free(alpha_one[i]);
		BN_clear_free(alpha_zero[i]);
		BN_clear_free(zi[i]);
		BN_clear_free(f_star[i]);
		for (size_t j=0; j<3; j++) {
			BN_clear_free(zij[i][j]);
			BN_clear_free(batch_y[i][j]);
			BN_clear_free(ytilde_ij[i][j]);
		}
	}
	for(size_t k=0; k<R; k++) {
		BN_clear_free(mu[k]);
		BN_clear_free(mutilde[k]);
		BN_clear_free(zeta[k]);
		BN_clear_free(tau[k]);
		BN_clear_free(dk[k]);
		BN_clear_free(fk[k]);
	}

	BN_CTX_free(bn_ctx);
	
	free(char_Dx);
	free(char_Dy);
	free(char_Dstar);
	for (size_t k=0; k<R; k++) {
		free(char_dk[k]);
	}
	free(char_Fx);
	free(char_Fy);
	free(char_Fstar);
	for (size_t k=0; k<R; k++) {
		free(char_fk[k]);
	}

	return 0;
}

// Generates a random point in the group
// Used to get random generators
void EC_POINT_generate_random_point(EC_GROUP *group, EC_POINT *p) {
	BN_CTX *bn_ctx = BN_CTX_new();
	BIGNUM *k = BN_new();

	EC_GROUP_get_order(group, k, bn_ctx);
	BN_rand(k, BN_num_bits(k), BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
	EC_POINT_mul(group, p, k, NULL, NULL, bn_ctx);

	BN_clear_free(k);
	BN_CTX_free(bn_ctx);
}


/********************************/ 
/*								*/
/* 			MAIN TEST	 		*/
/*								*/
/********************************/

int main(int argc, char **argv) {
	
	// init EC group G and generate 3 set of generators G0, Gi, Gi,j
	EC_GROUP *group_Gcom = EC_GROUP_new_by_curve_name(NID_secp256k1);
	
	EC_POINT *G_i[N], *G_ij[N][3];
	for (size_t i=0; i<N; i++) {
		G_i[i] = EC_POINT_new(group_Gcom);
		EC_POINT_generate_random_point(group_Gcom, G_i[i]);
		for (size_t j=0; j<3; j++) {
			G_ij[i][j] = EC_POINT_new(group_Gcom);
			EC_POINT_generate_random_point(group_Gcom, G_ij[i][j]);
		}
	}
	// init EC group Gtilte
	EC_POINT *Gtilde[R];
	for (size_t k=0; k<R; k++) {
		Gtilde[k] = EC_POINT_new(group_Gcom);
		EC_POINT_generate_random_point(group_Gcom, Gtilde[k]);
	}
	// init EC group H and generate generators Hi
	EC_GROUP *group_Hcom = EC_GROUP_new_by_curve_name(NID_secp256k1);
	EC_POINT *H_i[N];
	for (size_t i=0; i<N; i++) {
		H_i[i] = EC_POINT_new(group_Hcom);
		EC_POINT_generate_random_point(group_Hcom, H_i[i]);
	}

    // init integers x
    BIGNUM *x[N];
    for (size_t i=0; i<N; i++) {
        x[i] = BN_new();
        if (!BN_rand(x[i], B, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY))
            printf("Not enough entropy\n");
    }

	// init commitment x
	BN_CTX *bn_ctx = BN_CTX_new();
	EC_POINT *commitment_x = EC_POINT_new(group_Gcom);
	BIGNUM *r_x = BN_new();
	BN_rand(r_x, S, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
	EC_POINTs_mul(group_Gcom, commitment_x, r_x, N, (const EC_POINT **) G_i, (const BIGNUM **) x, bn_ctx);

	// Compute challenges of verifier (offline)	
	BIGNUM *scalar_gamma = BN_new();
	BN_rand(scalar_gamma, Omega, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);

	BIGNUM *batch_gamma[4][R][N]; 
	for (size_t j=0; j<4; j++) {
		for(size_t k=0; k<R; k++) {
			for (size_t i=0; i<N; i++) {
				batch_gamma[j][k][i] = BN_new();
				BN_rand(batch_gamma[j][k][i], Gamma, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
			}
		}
	}

	double timing[2];
	timing[0] = (double) 0;
	timing[1] = (double) 0;

	int nb_test = 1000;
	int nb_fail = 0;
	int res =0;

	// Run the protocol nb_test times and logs the number of aborts 
	// due to masking fails (aka rejection sampling).
	for (size_t i=0; i<nb_test; i++) {
	    res = protocol_poso(group_Gcom, group_Hcom, commitment_x, x, r_x, scalar_gamma, batch_gamma, G_i, G_ij, Gtilde, H_i, timing);
		if (res > 0) {
			nb_fail += 1;
		}
	}

	printf("Total number of rejected protocol: %d / %d\n", nb_fail, nb_test);

	printf("Total time taken by CPU for PROVE: %fs\n", timing[0]/(nb_test-nb_fail));

	printf("Total time taken by CPU for VERIF: %fs\n", timing[1]/(nb_test-nb_fail));

		
	for (size_t i=0; i<N; i++) {
		BN_clear_free(x[i]);
	}
	for (size_t j=0; j<4; j++) {
		for(size_t k=0; k<R; k++) {
			for (size_t i=0; i<N; i++) {
				BN_clear_free(batch_gamma[j][k][i]);
			}
		}
	}
	BN_clear_free(r_x);
	BN_clear_free(scalar_gamma);
	
	for (size_t i=0; i<N; i++) {
		EC_POINT_free(G_i[i]);
		EC_POINT_free(H_i[i]);
		for (size_t j=0; j<3; j++) {
			EC_POINT_free(G_ij[i][j]);
		}
	}
	for(size_t k=0; k<R; k++) {
		EC_POINT_free(Gtilde[k]);
	}

	EC_POINT_free(commitment_x);

	EC_GROUP_clear_free(group_Gcom);
	EC_GROUP_clear_free(group_Hcom);

	BN_CTX_free(bn_ctx);

    return 0;
}