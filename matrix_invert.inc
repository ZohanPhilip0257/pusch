#ifndef _MATRIX_INVERT_INC_
#define _MATRIX_INVERT_INC_


// Inverse = conj(input)/(abs(input)^2)
inline void InverseMatrix1x1HTranspose(VECT_T src_re[1][1], VECT_T src_im[1][1],
                                VECT_T dst_re[1][1], VECT_T dst_im[1][1])
{
    VECT_T invAbsSqr = VECT_MUL(src_re[0][0], src_re[0][0]);

#ifdef USE_FMA
    invAbsSqr = VECT_MUL_ADD(src_im[0][0], src_im[0][0], invAbsSqr);
#else
    VECT_T temp;
    temp = VECT_MUL(src_im[0][0], src_im[0][0]);
    invAbsSqr = VECT_ADD(invAbsSqr, temp);
#endif /* USE_FMA */

    invAbsSqr = VECT_RCP(invAbsSqr);

    dst_re[0][0] = VECT_MUL(src_re[0][0], invAbsSqr);
    dst_im[0][0] = VECT_MUL(VECT_NEG(src_im[0][0]), invAbsSqr);
}

/*void InverseMatrix2x2(VECT_T src_re[2][2], VECT_T src_im[2][2],
                      VECT_T dst_re[2][2], VECT_T dst_im[2][2])
{
    VECT_T ad_re;
    VECT_T ad_im;
    VECT_T bc_re;
    VECT_T bc_im;
    VECT_T det_re;
    VECT_T det_im;
    VECT_T temp;
    VECT_T temp2;
    int ii;
    int ij;

    // Create new matrix
    dst_re[0][0] = src_re[1][1];
    dst_im[0][0] = src_im[1][1];
    dst_re[1][1] = src_re[0][0];
    dst_im[1][1] = src_im[0][0];

    dst_re[0][1] = VECT_NEG(src_re[0][1]);
    dst_im[0][1] = VECT_NEG(src_im[0][1]);
    dst_re[1][0] = VECT_NEG(src_re[1][0]);
    dst_im[1][0] = VECT_NEG(src_im[1][0]);

    // Get determinant
#ifdef USE_FMA
    ad_re = VECT_MUL(src_re[0][0], src_re[1][1]);
    ad_re = VECT_NMUL_ADD(src_im[0][0], src_im[1][1] ,ad_re);
    ad_im = VECT_MUL(src_re[0][0], src_im[1][1]);
    ad_im = VECT_MUL_ADD(src_im[0][0], src_re[1][1], ad_im);

    bc_re = VECT_MUL(src_re[0][1], src_re[1][0]);
    bc_re = VECT_NMUL_ADD(src_im[0][1], src_im[1][0], bc_re);
    bc_im = VECT_MUL(src_re[0][1], src_im[1][0]);;
    bc_im = VECT_MUL_ADD(src_im[0][1], src_re[1][0], bc_im);
#else
    ad_re = VECT_MUL(src_re[0][0], src_re[1][1]);
    temp = VECT_MUL(src_im[0][0], src_im[1][1]);
    ad_re = VECT_SUB(ad_re, temp);
    ad_im = VECT_MUL(src_re[0][0], src_im[1][1]);
    temp = VECT_MUL(src_im[0][0], src_re[1][1]);
    ad_im = VECT_ADD(ad_im, temp);

    bc_re = VECT_MUL(src_re[0][1], src_re[1][0]);
    temp = VECT_MUL(src_im[0][1], src_im[1][0]);
    bc_re = VECT_SUB(bc_re, temp);
    bc_im = VECT_MUL(src_re[0][1], src_im[1][0]);
    temp = VECT_MUL(src_im[0][1], src_re[1][0]);
    bc_im = VECT_ADD(bc_im, temp);
#endif

    // Get the inverse
    det_re = VECT_SUB(ad_re, bc_re);
    det_im = VECT_SUB(ad_im, bc_im);

#ifdef USE_FMA
    temp = VECT_MUL(det_re, det_re);
    temp = VECT_MUL_ADD(det_im, det_im, temp);
#else
    temp = VECT_MUL(det_re, det_re);
    temp2 = VECT_MUL(det_im, det_im);
    temp = VECT_ADD(temp, temp2);
#endif

    temp = VECT_RCP(temp);
    det_re = VECT_MUL(det_re, temp);
    det_im = VECT_MUL(det_im, temp);
    det_im = VECT_NEG(det_im);

    // multiply across all cells
    for (ii = 0; ii < 2; ii++) {
        for (ij = 0; ij < 2; ij++) {
#ifdef USE_FMA
            ad_re = VECT_MUL(det_re, dst_re[ii][ij]);
            ad_re = VECT_NMUL_ADD(det_im, dst_im[ii][ij], ad_re);

            ad_im = VECT_MUL(det_re, dst_im[ii][ij]);
            dst_im[ii][ij] = VECT_MUL_ADD(det_im, dst_re[ii][ij], ad_im);
            dst_re[ii][ij] = ad_re;
#else
            ad_re = VECT_MUL(det_re, dst_re[ii][ij]);
            temp = VECT_MUL(det_im, dst_im[ii][ij]);
            ad_re = VECT_SUB(ad_re, temp);

            ad_im = VECT_MUL(det_re, dst_im[ii][ij]);
            temp = VECT_MUL(det_im, dst_re[ii][ij]);
            dst_im[ii][ij] = VECT_ADD(ad_im, temp);
            dst_re[ii][ij] = ad_re;
#endif
        }
    }
}
*/

inline void InverseMatrix2x2H(VECT_T src_re[2][2], VECT_T src_im[2][2],
                       VECT_T dst_re[2][2], VECT_T dst_im[2][2])
{
    VECT_T ad_re;
    // VECT_T ad_im;
    VECT_T bc_re;
    // VECT_T det_re;
    // VECT_T det_im;
    VECT_T temp;
    // VECT_T temp2;
    // int ii;
    // int ij;

    // Create new matrix
    dst_re[0][0] = src_re[1][1];
    dst_im[0][0] = src_im[1][1];
    dst_re[1][1] = src_re[0][0];
    dst_im[1][1] = src_im[0][0];

    dst_re[0][1] = VECT_NEG(src_re[0][1]);
    dst_im[0][1] = VECT_NEG(src_im[0][1]);

    // Get determinant
    ad_re = VECT_MUL(src_re[0][0], src_re[1][1]);
    bc_re = VECT_MUL(src_re[0][1], src_re[0][1]);
#ifdef USE_FMA
    bc_re = VECT_MUL_ADD(src_im[0][1], src_im[0][1], bc_re);
#else
    temp = VECT_MUL(src_im[0][1], src_im[0][1]);
    bc_re = VECT_ADD(bc_re, temp);
#endif /* USE_FMA */
    temp = VECT_SUB(ad_re, bc_re);

    // Get the inverse
    bc_re = VECT_RCP(temp);

    // multiply across all cells
    dst_re[0][0] = VECT_MUL(dst_re[0][0], bc_re);
    dst_re[0][1] = VECT_MUL(dst_re[0][1], bc_re);
    dst_im[0][1] = VECT_MUL(dst_im[0][1], bc_re);
    dst_re[1][1] = VECT_MUL(dst_re[1][1], bc_re);

    // Fill lower triangle
    dst_re[1][0] = dst_re[0][1];
    dst_im[1][0] = VECT_NEG(dst_im[0][1]);
}

// Stores result in transposed form
inline void InverseMatrix2x2HTranspose(VECT_T src_re[2][2], VECT_T src_im[2][2],
                                VECT_T dst_re[2][2], VECT_T dst_im[2][2])
{
    VECT_T ad_re;
    // VECT_T ad_im;
    VECT_T bc_re;
    // VECT_T det_re;
    // VECT_T det_im;
    VECT_T temp;
    // VECT_T temp2;
    // int ii;
    // int ij;

    // Create new matrix
    dst_re[0][0] = src_re[1][1];
    dst_im[0][0] = src_im[1][1];
    dst_re[1][1] = src_re[0][0];
    dst_im[1][1] = src_im[0][0];

    dst_re[0][1] = VECT_NEG(src_re[0][1]);
    dst_im[0][1] = VECT_NEG(src_im[0][1]);

    // Get determinant
    ad_re = VECT_MUL(src_re[0][0], src_re[1][1]);
    bc_re = VECT_MUL(src_re[0][1], src_re[0][1]);
#ifdef USE_FMA
    bc_re = VECT_MUL_ADD(src_im[0][1], src_im[0][1], bc_re);
#else
    temp = VECT_MUL(src_im[0][1], src_im[0][1]);
    bc_re = VECT_ADD(bc_re, temp);
#endif /* USE_FMA */
    temp = VECT_SUB(ad_re, bc_re);

    // Get the inverse
    bc_re = VECT_RCP(temp);

    // multiply across all cells
    dst_re[0][0] = VECT_MUL(dst_re[0][0], bc_re);
    dst_re[1][0] = VECT_MUL(dst_re[0][1], bc_re);
    dst_im[1][0] = VECT_MUL(dst_im[0][1], bc_re);
    dst_re[1][1] = VECT_MUL(dst_re[1][1], bc_re);

    // Fill lower triangle
    dst_re[0][1] = dst_re[1][0];
    dst_im[0][1] = VECT_NEG(dst_im[1][0]);
}

/*void InverseMatrix4x4(VECT_T src_re[4][4], VECT_T src_im[4][4],
                      VECT_T dst_re[4][4], VECT_T dst_im[4][4])
{
    VECT_T ma_re[2][2];
    VECT_T ma_im[2][2];
    VECT_T mb_re[2][2];
    VECT_T mb_im[2][2];
    VECT_T mc_re[2][2];
    VECT_T mc_im[2][2];
    VECT_T md_re[2][2];
    VECT_T md_im[2][2];

    VECT_T ma_out_re[2][2];
    VECT_T ma_out_im[2][2];
    VECT_T mb_out_re[2][2];
    VECT_T mb_out_im[2][2];
    VECT_T mc_out_re[2][2];
    VECT_T mc_out_im[2][2];
    VECT_T md_out_re[2][2];
    VECT_T md_out_im[2][2];

    VECT_T ma_inv_re[2][2];
    VECT_T ma_inv_im[2][2];

    VECT_T ma_temp_re[2][2];
    VECT_T ma_temp_im[2][2];

    VECT_T ma_temp2_re[2][2];
    VECT_T ma_temp2_im[2][2];

    VECT_T ma_cab_re[2][2];
    VECT_T ma_cab_im[2][2];

    VECT_T ma_ab_re[2][2];
    VECT_T ma_ab_im[2][2];

    VECT_T ma_ca_re[2][2];
    VECT_T ma_ca_im[2][2];

    int ii;
    int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];

    mb_re[0][0] = src_re[0][2]; mb_im[0][0] = src_im[0][2];
    mb_re[0][1] = src_re[0][3]; mb_im[0][1] = src_im[0][3];
    mb_re[1][0] = src_re[1][2]; mb_im[1][0] = src_im[1][2];
    mb_re[1][1] = src_re[1][3]; mb_im[1][1] = src_im[1][3];

    mc_re[0][0] = src_re[2][0]; mc_im[0][0] = src_im[2][0];
    mc_re[0][1] = src_re[2][1]; mc_im[0][1] = src_im[2][1];
    mc_re[1][0] = src_re[3][0]; mc_im[1][0] = src_im[3][0];
    mc_re[1][1] = src_re[3][1]; mc_im[1][1] = src_im[3][1];

    md_re[0][0] = src_re[2][2]; md_im[0][0] = src_im[2][2];
    md_re[0][1] = src_re[2][3]; md_im[0][1] = src_im[2][3];
    md_re[1][0] = src_re[3][2]; md_im[1][0] = src_im[3][2];
    md_re[1][1] = src_re[3][3]; md_im[1][1] = src_im[3][3];

    // Get inv(A)
    InverseMatrix2x2(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMUL(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
              2, 2, 2);

    // C*inv(A)*B
    MATRIXMUL(ma_ca_re, ma_ca_im, mb_re, mb_im, ma_temp2_re, ma_temp2_im,
              2, 2, 2);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              2, 2);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix2x2(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMUL(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
              2, 2, 2);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMUL(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
              ma_temp_re, ma_temp_im, 2, 2, 2);
    // [B]
    MATRIXSUB(g_mzero2x2, g_mzero2x2, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 2, 2);

    // (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(md_out_re, md_out_im, ma_ca_re, ma_ca_im,
              ma_temp_re, ma_temp_im, 2, 2, 2);

    // [C]
    MATRIXSUB(g_mzero2x2, g_mzero2x2, ma_temp_re, ma_temp_im,
              mc_out_re, mc_out_im, 2, 2);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(ma_ab_re, ma_ab_im, ma_temp_re, ma_temp_im,
              ma_temp2_re, ma_temp2_im, 2, 2, 2);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 2, 2);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[0][1] = ma_out_re[0][1]; dst_im[0][1] = ma_out_im[0][1];
    dst_re[1][0] = ma_out_re[1][0]; dst_im[1][0] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];

    dst_re[0][2] = mb_out_re[0][0]; dst_im[0][2] = mb_out_im[0][0];
    dst_re[0][3] = mb_out_re[0][1]; dst_im[0][3] = mb_out_im[0][1];
    dst_re[1][2] = mb_out_re[1][0]; dst_im[1][2] = mb_out_im[1][0];
    dst_re[1][3] = mb_out_re[1][1]; dst_im[1][3] = mb_out_im[1][1];

    dst_re[2][0] = mc_out_re[0][0]; dst_im[2][0] = mc_out_im[0][0];
    dst_re[2][1] = mc_out_re[0][1]; dst_im[2][1] = mc_out_im[0][1];
    dst_re[3][0] = mc_out_re[1][0]; dst_im[3][0] = mc_out_im[1][0];
    dst_re[3][1] = mc_out_re[1][1]; dst_im[3][1] = mc_out_im[1][1];

    dst_re[2][2] = md_out_re[0][0]; dst_im[2][2] = md_out_im[0][0];
    dst_re[2][3] = md_out_re[0][1]; dst_im[2][3] = md_out_im[0][1];
    dst_re[3][2] = md_out_re[1][0]; dst_im[3][2] = md_out_im[1][0];
    dst_re[3][3] = md_out_re[1][1]; dst_im[3][3] = md_out_im[1][1];
}
*/

inline void InverseMatrix4x4H(VECT_T src_re[4][4], VECT_T src_im[4][4],VECT_T dst_re[4][4], VECT_T dst_im[4][4])
{
    VECT_T ma_re[2][2];
    VECT_T ma_im[2][2];
    VECT_T mb_re[2][2];
    VECT_T mb_im[2][2];
    VECT_T mc_re[2][2];
    VECT_T mc_im[2][2];
    VECT_T md_re[2][2];
    VECT_T md_im[2][2];

    VECT_T ma_out_re[2][2];
    VECT_T ma_out_im[2][2];
    VECT_T mb_out_re[2][2];
    VECT_T mb_out_im[2][2];
    VECT_T md_out_re[2][2];
    VECT_T md_out_im[2][2];

    VECT_T ma_inv_re[2][2];
    VECT_T ma_inv_im[2][2];

    VECT_T ma_temp_re[2][2];
    VECT_T ma_temp_im[2][2];

    VECT_T ma_temp2_re[2][2];
    VECT_T ma_temp2_im[2][2];

    VECT_T ma_cab_re[2][2];
    VECT_T ma_cab_im[2][2];

    VECT_T ma_ab_re[2][2];
    VECT_T ma_ab_im[2][2];

    VECT_T ma_ca_re[2][2];
    VECT_T ma_ca_im[2][2];

    // int ii;
    // int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];

    mb_re[0][0] = src_re[0][2]; mb_im[0][0] = src_im[0][2];
    mb_re[0][1] = src_re[0][3]; mb_im[0][1] = src_im[0][3];
    mb_re[1][0] = src_re[1][2]; mb_im[1][0] = src_im[1][2];
    mb_re[1][1] = src_re[1][3]; mb_im[1][1] = src_im[1][3];

    mc_re[0][0] = src_re[2][0]; mc_im[0][0] = src_im[2][0];
    mc_re[0][1] = src_re[2][1]; mc_im[0][1] = src_im[2][1];
    mc_re[1][0] = src_re[3][0]; mc_im[1][0] = src_im[3][0];
    mc_re[1][1] = src_re[3][1]; mc_im[1][1] = src_im[3][1];

    md_re[0][0] = src_re[2][2]; md_im[0][0] = src_im[2][2];
    md_re[0][1] = src_re[2][3]; md_im[0][1] = src_im[2][3];
    md_re[1][0] = src_re[3][2]; md_im[1][0] = src_im[3][2];
    md_re[1][1] = src_re[3][3]; md_im[1][1] = src_im[3][3];

    // Get inv(A)
    InverseMatrix2x2H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                2, 2, 2);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 2, 2, 2);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              2, 2);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix2x2H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                2, 2, 2);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 2, 2, 2);
    // [B]
    MATRIXSUB(g_mzero2x2, g_mzero2x2, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 2, 2);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 2, 2, 2);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 2, 2);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[0][1] = ma_out_re[0][1]; dst_im[0][1] = ma_out_im[0][1];
    dst_re[1][0] = ma_out_re[1][0]; dst_im[1][0] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];

    dst_re[0][2] = mb_out_re[0][0]; dst_im[0][2] = mb_out_im[0][0];
    dst_re[0][3] = mb_out_re[0][1]; dst_im[0][3] = mb_out_im[0][1];
    dst_re[1][2] = mb_out_re[1][0]; dst_im[1][2] = mb_out_im[1][0];
    dst_re[1][3] = mb_out_re[1][1]; dst_im[1][3] = mb_out_im[1][1];

    dst_re[2][0] = mb_out_re[0][0]; dst_im[2][0] = VECT_NEG(mb_out_im[0][0]);
    dst_re[2][1] = mb_out_re[1][0]; dst_im[2][1] = VECT_NEG(mb_out_im[1][0]);
    dst_re[3][0] = mb_out_re[0][1]; dst_im[3][0] = VECT_NEG(mb_out_im[0][1]);
    dst_re[3][1] = mb_out_re[1][1]; dst_im[3][1] = VECT_NEG(mb_out_im[1][1]);

    dst_re[2][2] = md_out_re[0][0]; dst_im[2][2] = md_out_im[0][0];
    dst_re[2][3] = md_out_re[0][1]; dst_im[2][3] = md_out_im[0][1];
    dst_re[3][2] = md_out_re[1][0]; dst_im[3][2] = md_out_im[1][0];
    dst_re[3][3] = md_out_re[1][1]; dst_im[3][3] = md_out_im[1][1];
}

// Stores result in transposed form
inline void InverseMatrix4x4HTranspose(VECT_T src_re[4][4], VECT_T src_im[4][4],
                                VECT_T dst_re[4][4], VECT_T dst_im[4][4])
{
    VECT_T ma_re[2][2];
    VECT_T ma_im[2][2];
    VECT_T mb_re[2][2];
    VECT_T mb_im[2][2];
    VECT_T mc_re[2][2];
    VECT_T mc_im[2][2];
    VECT_T md_re[2][2];
    VECT_T md_im[2][2];

    VECT_T ma_out_re[2][2];
    VECT_T ma_out_im[2][2];
    VECT_T mb_out_re[2][2];
    VECT_T mb_out_im[2][2];
    VECT_T md_out_re[2][2];
    VECT_T md_out_im[2][2];

    VECT_T ma_inv_re[2][2];
    VECT_T ma_inv_im[2][2];

    VECT_T ma_temp_re[2][2];
    VECT_T ma_temp_im[2][2];

    VECT_T ma_temp2_re[2][2];
    VECT_T ma_temp2_im[2][2];

    VECT_T ma_cab_re[2][2];
    VECT_T ma_cab_im[2][2];

    VECT_T ma_ab_re[2][2];
    VECT_T ma_ab_im[2][2];

    VECT_T ma_ca_re[2][2];
    VECT_T ma_ca_im[2][2];

    // int ii;
    // int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];

    mb_re[0][0] = src_re[0][2]; mb_im[0][0] = src_im[0][2];
    mb_re[0][1] = src_re[0][3]; mb_im[0][1] = src_im[0][3];
    mb_re[1][0] = src_re[1][2]; mb_im[1][0] = src_im[1][2];
    mb_re[1][1] = src_re[1][3]; mb_im[1][1] = src_im[1][3];

    mc_re[0][0] = src_re[2][0]; mc_im[0][0] = src_im[2][0];
    mc_re[0][1] = src_re[2][1]; mc_im[0][1] = src_im[2][1];
    mc_re[1][0] = src_re[3][0]; mc_im[1][0] = src_im[3][0];
    mc_re[1][1] = src_re[3][1]; mc_im[1][1] = src_im[3][1];

    md_re[0][0] = src_re[2][2]; md_im[0][0] = src_im[2][2];
    md_re[0][1] = src_re[2][3]; md_im[0][1] = src_im[2][3];
    md_re[1][0] = src_re[3][2]; md_im[1][0] = src_im[3][2];
    md_re[1][1] = src_re[3][3]; md_im[1][1] = src_im[3][3];

    // Get inv(A)
    InverseMatrix2x2H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                2, 2, 2);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 2, 2, 2);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              2, 2);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix2x2H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                2, 2, 2);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 2, 2, 2);
    // [B]
    MATRIXSUB(g_mzero2x2, g_mzero2x2, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 2, 2);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 2, 2, 2);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 2, 2);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[1][0] = ma_out_re[0][1]; dst_im[1][0] = ma_out_im[0][1];
    dst_re[0][1] = ma_out_re[1][0]; dst_im[0][1] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];

    dst_re[2][0] = mb_out_re[0][0]; dst_im[2][0] = mb_out_im[0][0];
    dst_re[3][0] = mb_out_re[0][1]; dst_im[3][0] = mb_out_im[0][1];
    dst_re[2][1] = mb_out_re[1][0]; dst_im[2][1] = mb_out_im[1][0];
    dst_re[3][1] = mb_out_re[1][1]; dst_im[3][1] = mb_out_im[1][1];

    dst_re[0][2] = mb_out_re[0][0]; dst_im[0][2] = VECT_NEG(mb_out_im[0][0]);
    dst_re[1][2] = mb_out_re[1][0]; dst_im[1][2] = VECT_NEG(mb_out_im[1][0]);
    dst_re[0][3] = mb_out_re[0][1]; dst_im[0][3] = VECT_NEG(mb_out_im[0][1]);
    dst_re[1][3] = mb_out_re[1][1]; dst_im[1][3] = VECT_NEG(mb_out_im[1][1]);

    dst_re[2][2] = md_out_re[0][0]; dst_im[2][2] = md_out_im[0][0];
    dst_re[3][2] = md_out_re[0][1]; dst_im[3][2] = md_out_im[0][1];
    dst_re[2][3] = md_out_re[1][0]; dst_im[2][3] = md_out_im[1][0];
    dst_re[3][3] = md_out_re[1][1]; dst_im[3][3] = md_out_im[1][1];
}

/*void InverseMatrix8x8(VECT_T src_re[8][8], VECT_T src_im[8][8],
                      VECT_T dst_re[8][8], VECT_T dst_im[8][8])
{
    VECT_T ma_re[4][4];
    VECT_T ma_im[4][4];
    VECT_T mb_re[4][4];
    VECT_T mb_im[4][4];
    VECT_T mc_re[4][4];
    VECT_T mc_im[4][4];
    VECT_T md_re[4][4];
    VECT_T md_im[4][4];

    VECT_T ma_out_re[4][4];
    VECT_T ma_out_im[4][4];
    VECT_T mb_out_re[4][4];
    VECT_T mb_out_im[4][4];
    VECT_T mc_out_re[4][4];
    VECT_T mc_out_im[4][4];
    VECT_T md_out_re[4][4];
    VECT_T md_out_im[4][4];

    VECT_T ma_inv_re[4][4];
    VECT_T ma_inv_im[4][4];

    VECT_T ma_temp_re[4][4];
    VECT_T ma_temp_im[4][4];

    VECT_T ma_temp2_re[4][4];
    VECT_T ma_temp2_im[4][4];

    VECT_T ma_cab_re[4][4];
    VECT_T ma_cab_im[4][4];

    VECT_T ma_ab_re[4][4];
    VECT_T ma_ab_im[4][4];

    VECT_T ma_ca_re[4][4];
    VECT_T ma_ca_im[4][4];

    int ii;
    int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[0][2] = src_re[0][2]; ma_im[0][2] = src_im[0][2];
    ma_re[0][3] = src_re[0][3]; ma_im[0][3] = src_im[0][3];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];
    ma_re[1][2] = src_re[1][2]; ma_im[1][2] = src_im[1][2];
    ma_re[1][3] = src_re[1][3]; ma_im[1][3] = src_im[1][3];
    ma_re[2][0] = src_re[2][0]; ma_im[2][0] = src_im[2][0];
    ma_re[2][1] = src_re[2][1]; ma_im[2][1] = src_im[2][1];
    ma_re[2][2] = src_re[2][2]; ma_im[2][2] = src_im[2][2];
    ma_re[2][3] = src_re[2][3]; ma_im[2][3] = src_im[2][3];
    ma_re[3][0] = src_re[3][0]; ma_im[3][0] = src_im[3][0];
    ma_re[3][1] = src_re[3][1]; ma_im[3][1] = src_im[3][1];
    ma_re[3][2] = src_re[3][2]; ma_im[3][2] = src_im[3][2];
    ma_re[3][3] = src_re[3][3]; ma_im[3][3] = src_im[3][3];

    mb_re[0][0] = src_re[0][4]; mb_im[0][0] = src_im[0][4];
    mb_re[0][1] = src_re[0][5]; mb_im[0][1] = src_im[0][5];
    mb_re[0][2] = src_re[0][6]; mb_im[0][2] = src_im[0][6];
    mb_re[0][3] = src_re[0][7]; mb_im[0][3] = src_im[0][7];
    mb_re[1][0] = src_re[1][4]; mb_im[1][0] = src_im[1][4];
    mb_re[1][1] = src_re[1][5]; mb_im[1][1] = src_im[1][5];
    mb_re[1][2] = src_re[1][6]; mb_im[1][2] = src_im[1][6];
    mb_re[1][3] = src_re[1][7]; mb_im[1][3] = src_im[1][7];
    mb_re[2][0] = src_re[2][4]; mb_im[2][0] = src_im[2][4];
    mb_re[2][1] = src_re[2][5]; mb_im[2][1] = src_im[2][5];
    mb_re[2][2] = src_re[2][6]; mb_im[2][2] = src_im[2][6];
    mb_re[2][3] = src_re[2][7]; mb_im[2][3] = src_im[2][7];
    mb_re[3][0] = src_re[3][4]; mb_im[3][0] = src_im[3][4];
    mb_re[3][1] = src_re[3][5]; mb_im[3][1] = src_im[3][5];
    mb_re[3][2] = src_re[3][6]; mb_im[3][2] = src_im[3][6];
    mb_re[3][3] = src_re[3][7]; mb_im[3][3] = src_im[3][7];

    mc_re[0][0] = src_re[4][0]; mc_im[0][0] = src_im[4][0];
    mc_re[0][1] = src_re[4][1]; mc_im[0][1] = src_im[4][1];
    mc_re[0][2] = src_re[4][2]; mc_im[0][2] = src_im[4][2];
    mc_re[0][3] = src_re[4][3]; mc_im[0][3] = src_im[4][3];
    mc_re[1][0] = src_re[5][0]; mc_im[1][0] = src_im[5][0];
    mc_re[1][1] = src_re[5][1]; mc_im[1][1] = src_im[5][1];
    mc_re[1][2] = src_re[5][2]; mc_im[1][2] = src_im[5][2];
    mc_re[1][3] = src_re[5][3]; mc_im[1][3] = src_im[5][3];
    mc_re[2][0] = src_re[6][0]; mc_im[2][0] = src_im[6][0];
    mc_re[2][1] = src_re[6][1]; mc_im[2][1] = src_im[6][1];
    mc_re[2][2] = src_re[6][2]; mc_im[2][2] = src_im[6][2];
    mc_re[2][3] = src_re[6][3]; mc_im[2][3] = src_im[6][3];
    mc_re[3][0] = src_re[7][0]; mc_im[3][0] = src_im[7][0];
    mc_re[3][1] = src_re[7][1]; mc_im[3][1] = src_im[7][1];
    mc_re[3][2] = src_re[7][2]; mc_im[3][2] = src_im[7][2];
    mc_re[3][3] = src_re[7][3]; mc_im[3][3] = src_im[7][3];

    md_re[0][0] = src_re[4][4]; md_im[0][0] = src_im[4][4];
    md_re[0][1] = src_re[4][5]; md_im[0][1] = src_im[4][5];
    md_re[0][2] = src_re[4][6]; md_im[0][2] = src_im[4][6];
    md_re[0][3] = src_re[4][7]; md_im[0][3] = src_im[4][7];
    md_re[1][0] = src_re[5][4]; md_im[1][0] = src_im[5][4];
    md_re[1][1] = src_re[5][5]; md_im[1][1] = src_im[5][5];
    md_re[1][2] = src_re[5][6]; md_im[1][2] = src_im[5][6];
    md_re[1][3] = src_re[5][7]; md_im[1][3] = src_im[5][7];
    md_re[2][0] = src_re[6][4]; md_im[2][0] = src_im[6][4];
    md_re[2][1] = src_re[6][5]; md_im[2][1] = src_im[6][5];
    md_re[2][2] = src_re[6][6]; md_im[2][2] = src_im[6][6];
    md_re[2][3] = src_re[6][7]; md_im[2][3] = src_im[6][7];
    md_re[3][0] = src_re[7][4]; md_im[3][0] = src_im[7][4];
    md_re[3][1] = src_re[7][5]; md_im[3][1] = src_im[7][5];
    md_re[3][2] = src_re[7][6]; md_im[3][2] = src_im[7][6];
    md_re[3][3] = src_re[7][7]; md_im[3][3] = src_im[7][7];

    // Get inv(A)
    InverseMatrix4x4(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMUL(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
              4, 4, 4);

    // C*inv(A)*B
    MATRIXMUL(ma_ca_re, ma_ca_im, mb_re, mb_im, ma_temp2_re, ma_temp2_im,
              4, 4, 4);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              4, 4);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix4x4(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMUL(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
              4, 4, 4);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMUL(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
              ma_temp_re, ma_temp_im, 4, 4, 4);
    // [B]
    MATRIXSUB(g_mzero4x4, g_mzero4x4, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 4, 4);

    // (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(md_out_re, md_out_im, ma_ca_re, ma_ca_im,
              ma_temp_re, ma_temp_im, 4, 4, 4);

    // [C]
    MATRIXSUB(g_mzero4x4, g_mzero4x4, ma_temp_re, ma_temp_im,
              mc_out_re, mc_out_im, 4, 4);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(ma_ab_re, ma_ab_im, ma_temp_re, ma_temp_im,
              ma_temp2_re, ma_temp2_im, 4, 4, 4);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 4, 4);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[0][1] = ma_out_re[0][1]; dst_im[0][1] = ma_out_im[0][1];
    dst_re[0][2] = ma_out_re[0][2]; dst_im[0][2] = ma_out_im[0][2];
    dst_re[0][3] = ma_out_re[0][3]; dst_im[0][3] = ma_out_im[0][3];
    dst_re[1][0] = ma_out_re[1][0]; dst_im[1][0] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];
    dst_re[1][2] = ma_out_re[1][2]; dst_im[1][2] = ma_out_im[1][2];
    dst_re[1][3] = ma_out_re[1][3]; dst_im[1][3] = ma_out_im[1][3];
    dst_re[2][0] = ma_out_re[2][0]; dst_im[2][0] = ma_out_im[2][0];
    dst_re[2][1] = ma_out_re[2][1]; dst_im[2][1] = ma_out_im[2][1];
    dst_re[2][2] = ma_out_re[2][2]; dst_im[2][2] = ma_out_im[2][2];
    dst_re[2][3] = ma_out_re[2][3]; dst_im[2][3] = ma_out_im[2][3];
    dst_re[3][0] = ma_out_re[3][0]; dst_im[3][0] = ma_out_im[3][0];
    dst_re[3][1] = ma_out_re[3][1]; dst_im[3][1] = ma_out_im[3][1];
    dst_re[3][2] = ma_out_re[3][2]; dst_im[3][2] = ma_out_im[3][2];
    dst_re[3][3] = ma_out_re[3][3]; dst_im[3][3] = ma_out_im[3][3];

    dst_re[0][4] = mb_out_re[0][0]; dst_im[0][4] = mb_out_im[0][0];
    dst_re[0][5] = mb_out_re[0][1]; dst_im[0][5] = mb_out_im[0][1];
    dst_re[0][6] = mb_out_re[0][2]; dst_im[0][6] = mb_out_im[0][2];
    dst_re[0][7] = mb_out_re[0][3]; dst_im[0][7] = mb_out_im[0][3];
    dst_re[1][4] = mb_out_re[1][0]; dst_im[1][4] = mb_out_im[1][0];
    dst_re[1][5] = mb_out_re[1][1]; dst_im[1][5] = mb_out_im[1][1];
    dst_re[1][6] = mb_out_re[1][2]; dst_im[1][6] = mb_out_im[1][2];
    dst_re[1][7] = mb_out_re[1][3]; dst_im[1][7] = mb_out_im[1][3];
    dst_re[2][4] = mb_out_re[2][0]; dst_im[2][4] = mb_out_im[2][0];
    dst_re[2][5] = mb_out_re[2][1]; dst_im[2][5] = mb_out_im[2][1];
    dst_re[2][6] = mb_out_re[2][2]; dst_im[2][6] = mb_out_im[2][2];
    dst_re[2][7] = mb_out_re[2][3]; dst_im[2][7] = mb_out_im[2][3];
    dst_re[3][4] = mb_out_re[3][0]; dst_im[3][4] = mb_out_im[3][0];
    dst_re[3][5] = mb_out_re[3][1]; dst_im[3][5] = mb_out_im[3][1];
    dst_re[3][6] = mb_out_re[3][2]; dst_im[3][6] = mb_out_im[3][2];
    dst_re[3][7] = mb_out_re[3][3]; dst_im[3][7] = mb_out_im[3][3];

    dst_re[4][0] = mc_out_re[0][0]; dst_im[4][0] = mc_out_im[0][0];
    dst_re[4][1] = mc_out_re[0][1]; dst_im[4][1] = mc_out_im[0][1];
    dst_re[4][2] = mc_out_re[0][2]; dst_im[4][2] = mc_out_im[0][2];
    dst_re[4][3] = mc_out_re[0][3]; dst_im[4][3] = mc_out_im[0][3];
    dst_re[5][0] = mc_out_re[1][0]; dst_im[5][0] = mc_out_im[1][0];
    dst_re[5][1] = mc_out_re[1][1]; dst_im[5][1] = mc_out_im[1][1];
    dst_re[5][2] = mc_out_re[1][2]; dst_im[5][2] = mc_out_im[1][2];
    dst_re[5][3] = mc_out_re[1][3]; dst_im[5][3] = mc_out_im[1][3];
    dst_re[6][0] = mc_out_re[2][0]; dst_im[6][0] = mc_out_im[2][0];
    dst_re[6][1] = mc_out_re[2][1]; dst_im[6][1] = mc_out_im[2][1];
    dst_re[6][2] = mc_out_re[2][2]; dst_im[6][2] = mc_out_im[2][2];
    dst_re[6][3] = mc_out_re[2][3]; dst_im[6][3] = mc_out_im[2][3];
    dst_re[7][0] = mc_out_re[3][0]; dst_im[7][0] = mc_out_im[3][0];
    dst_re[7][1] = mc_out_re[3][1]; dst_im[7][1] = mc_out_im[3][1];
    dst_re[7][2] = mc_out_re[3][2]; dst_im[7][2] = mc_out_im[3][2];
    dst_re[7][3] = mc_out_re[3][3]; dst_im[7][3] = mc_out_im[3][3];

    dst_re[4][4] = md_out_re[0][0]; dst_im[4][4] = md_out_im[0][0];
    dst_re[4][5] = md_out_re[0][1]; dst_im[4][5] = md_out_im[0][1];
    dst_re[4][6] = md_out_re[0][2]; dst_im[4][6] = md_out_im[0][2];
    dst_re[4][7] = md_out_re[0][3]; dst_im[4][7] = md_out_im[0][3];
    dst_re[5][4] = md_out_re[1][0]; dst_im[5][4] = md_out_im[1][0];
    dst_re[5][5] = md_out_re[1][1]; dst_im[5][5] = md_out_im[1][1];
    dst_re[5][6] = md_out_re[1][2]; dst_im[5][6] = md_out_im[1][2];
    dst_re[5][7] = md_out_re[1][3]; dst_im[5][7] = md_out_im[1][3];
    dst_re[6][4] = md_out_re[2][0]; dst_im[6][4] = md_out_im[2][0];
    dst_re[6][5] = md_out_re[2][1]; dst_im[6][5] = md_out_im[2][1];
    dst_re[6][6] = md_out_re[2][2]; dst_im[6][6] = md_out_im[2][2];
    dst_re[6][7] = md_out_re[2][3]; dst_im[6][7] = md_out_im[2][3];
    dst_re[7][4] = md_out_re[3][0]; dst_im[7][4] = md_out_im[3][0];
    dst_re[7][5] = md_out_re[3][1]; dst_im[7][5] = md_out_im[3][1];
    dst_re[7][6] = md_out_re[3][2]; dst_im[7][6] = md_out_im[3][2];
    dst_re[7][7] = md_out_re[3][3]; dst_im[7][7] = md_out_im[3][3];

}
*/

inline void InverseMatrix8x8H(VECT_T src_re[8][8], VECT_T src_im[8][8],
                       VECT_T dst_re[8][8], VECT_T dst_im[8][8])
{
    VECT_T ma_re[4][4];
    VECT_T ma_im[4][4];
    VECT_T mb_re[4][4];
    VECT_T mb_im[4][4];
    VECT_T mc_re[4][4];
    VECT_T mc_im[4][4];
    VECT_T md_re[4][4];
    VECT_T md_im[4][4];

    VECT_T ma_out_re[4][4];
    VECT_T ma_out_im[4][4];
    VECT_T mb_out_re[4][4];
    VECT_T mb_out_im[4][4];
    VECT_T md_out_re[4][4];
    VECT_T md_out_im[4][4];

    VECT_T ma_inv_re[4][4];
    VECT_T ma_inv_im[4][4];

    VECT_T ma_temp_re[4][4];
    VECT_T ma_temp_im[4][4];

    VECT_T ma_temp2_re[4][4];
    VECT_T ma_temp2_im[4][4];

    VECT_T ma_cab_re[4][4];
    VECT_T ma_cab_im[4][4];

    VECT_T ma_ab_re[4][4];
    VECT_T ma_ab_im[4][4];

    VECT_T ma_ca_re[4][4];
    VECT_T ma_ca_im[4][4];

    // int ii;
    // int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[0][2] = src_re[0][2]; ma_im[0][2] = src_im[0][2];
    ma_re[0][3] = src_re[0][3]; ma_im[0][3] = src_im[0][3];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];
    ma_re[1][2] = src_re[1][2]; ma_im[1][2] = src_im[1][2];
    ma_re[1][3] = src_re[1][3]; ma_im[1][3] = src_im[1][3];
    ma_re[2][0] = src_re[2][0]; ma_im[2][0] = src_im[2][0];
    ma_re[2][1] = src_re[2][1]; ma_im[2][1] = src_im[2][1];
    ma_re[2][2] = src_re[2][2]; ma_im[2][2] = src_im[2][2];
    ma_re[2][3] = src_re[2][3]; ma_im[2][3] = src_im[2][3];
    ma_re[3][0] = src_re[3][0]; ma_im[3][0] = src_im[3][0];
    ma_re[3][1] = src_re[3][1]; ma_im[3][1] = src_im[3][1];
    ma_re[3][2] = src_re[3][2]; ma_im[3][2] = src_im[3][2];
    ma_re[3][3] = src_re[3][3]; ma_im[3][3] = src_im[3][3];

    mb_re[0][0] = src_re[0][4]; mb_im[0][0] = src_im[0][4];
    mb_re[0][1] = src_re[0][5]; mb_im[0][1] = src_im[0][5];
    mb_re[0][2] = src_re[0][6]; mb_im[0][2] = src_im[0][6];
    mb_re[0][3] = src_re[0][7]; mb_im[0][3] = src_im[0][7];
    mb_re[1][0] = src_re[1][4]; mb_im[1][0] = src_im[1][4];
    mb_re[1][1] = src_re[1][5]; mb_im[1][1] = src_im[1][5];
    mb_re[1][2] = src_re[1][6]; mb_im[1][2] = src_im[1][6];
    mb_re[1][3] = src_re[1][7]; mb_im[1][3] = src_im[1][7];
    mb_re[2][0] = src_re[2][4]; mb_im[2][0] = src_im[2][4];
    mb_re[2][1] = src_re[2][5]; mb_im[2][1] = src_im[2][5];
    mb_re[2][2] = src_re[2][6]; mb_im[2][2] = src_im[2][6];
    mb_re[2][3] = src_re[2][7]; mb_im[2][3] = src_im[2][7];
    mb_re[3][0] = src_re[3][4]; mb_im[3][0] = src_im[3][4];
    mb_re[3][1] = src_re[3][5]; mb_im[3][1] = src_im[3][5];
    mb_re[3][2] = src_re[3][6]; mb_im[3][2] = src_im[3][6];
    mb_re[3][3] = src_re[3][7]; mb_im[3][3] = src_im[3][7];

    mc_re[0][0] = src_re[4][0]; mc_im[0][0] = src_im[4][0];
    mc_re[0][1] = src_re[4][1]; mc_im[0][1] = src_im[4][1];
    mc_re[0][2] = src_re[4][2]; mc_im[0][2] = src_im[4][2];
    mc_re[0][3] = src_re[4][3]; mc_im[0][3] = src_im[4][3];
    mc_re[1][0] = src_re[5][0]; mc_im[1][0] = src_im[5][0];
    mc_re[1][1] = src_re[5][1]; mc_im[1][1] = src_im[5][1];
    mc_re[1][2] = src_re[5][2]; mc_im[1][2] = src_im[5][2];
    mc_re[1][3] = src_re[5][3]; mc_im[1][3] = src_im[5][3];
    mc_re[2][0] = src_re[6][0]; mc_im[2][0] = src_im[6][0];
    mc_re[2][1] = src_re[6][1]; mc_im[2][1] = src_im[6][1];
    mc_re[2][2] = src_re[6][2]; mc_im[2][2] = src_im[6][2];
    mc_re[2][3] = src_re[6][3]; mc_im[2][3] = src_im[6][3];
    mc_re[3][0] = src_re[7][0]; mc_im[3][0] = src_im[7][0];
    mc_re[3][1] = src_re[7][1]; mc_im[3][1] = src_im[7][1];
    mc_re[3][2] = src_re[7][2]; mc_im[3][2] = src_im[7][2];
    mc_re[3][3] = src_re[7][3]; mc_im[3][3] = src_im[7][3];

    md_re[0][0] = src_re[4][4]; md_im[0][0] = src_im[4][4];
    md_re[0][1] = src_re[4][5]; md_im[0][1] = src_im[4][5];
    md_re[0][2] = src_re[4][6]; md_im[0][2] = src_im[4][6];
    md_re[0][3] = src_re[4][7]; md_im[0][3] = src_im[4][7];
    md_re[1][0] = src_re[5][4]; md_im[1][0] = src_im[5][4];
    md_re[1][1] = src_re[5][5]; md_im[1][1] = src_im[5][5];
    md_re[1][2] = src_re[5][6]; md_im[1][2] = src_im[5][6];
    md_re[1][3] = src_re[5][7]; md_im[1][3] = src_im[5][7];
    md_re[2][0] = src_re[6][4]; md_im[2][0] = src_im[6][4];
    md_re[2][1] = src_re[6][5]; md_im[2][1] = src_im[6][5];
    md_re[2][2] = src_re[6][6]; md_im[2][2] = src_im[6][6];
    md_re[2][3] = src_re[6][7]; md_im[2][3] = src_im[6][7];
    md_re[3][0] = src_re[7][4]; md_im[3][0] = src_im[7][4];
    md_re[3][1] = src_re[7][5]; md_im[3][1] = src_im[7][5];
    md_re[3][2] = src_re[7][6]; md_im[3][2] = src_im[7][6];
    md_re[3][3] = src_re[7][7]; md_im[3][3] = src_im[7][7];

    // Get inv(A)
    InverseMatrix4x4H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                4, 4, 4);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 4, 4, 4);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              4, 4);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix4x4H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                4, 4, 4);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 4, 4, 4);
    // [B]
    MATRIXSUB(g_mzero4x4, g_mzero4x4, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 4, 4);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 4, 4, 4);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 4, 4);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[0][1] = ma_out_re[0][1]; dst_im[0][1] = ma_out_im[0][1];
    dst_re[0][2] = ma_out_re[0][2]; dst_im[0][2] = ma_out_im[0][2];
    dst_re[0][3] = ma_out_re[0][3]; dst_im[0][3] = ma_out_im[0][3];
    dst_re[1][0] = ma_out_re[1][0]; dst_im[1][0] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];
    dst_re[1][2] = ma_out_re[1][2]; dst_im[1][2] = ma_out_im[1][2];
    dst_re[1][3] = ma_out_re[1][3]; dst_im[1][3] = ma_out_im[1][3];
    dst_re[2][0] = ma_out_re[2][0]; dst_im[2][0] = ma_out_im[2][0];
    dst_re[2][1] = ma_out_re[2][1]; dst_im[2][1] = ma_out_im[2][1];
    dst_re[2][2] = ma_out_re[2][2]; dst_im[2][2] = ma_out_im[2][2];
    dst_re[2][3] = ma_out_re[2][3]; dst_im[2][3] = ma_out_im[2][3];
    dst_re[3][0] = ma_out_re[3][0]; dst_im[3][0] = ma_out_im[3][0];
    dst_re[3][1] = ma_out_re[3][1]; dst_im[3][1] = ma_out_im[3][1];
    dst_re[3][2] = ma_out_re[3][2]; dst_im[3][2] = ma_out_im[3][2];
    dst_re[3][3] = ma_out_re[3][3]; dst_im[3][3] = ma_out_im[3][3];

    dst_re[0][4] = mb_out_re[0][0]; dst_im[0][4] = mb_out_im[0][0];
    dst_re[0][5] = mb_out_re[0][1]; dst_im[0][5] = mb_out_im[0][1];
    dst_re[0][6] = mb_out_re[0][2]; dst_im[0][6] = mb_out_im[0][2];
    dst_re[0][7] = mb_out_re[0][3]; dst_im[0][7] = mb_out_im[0][3];
    dst_re[1][4] = mb_out_re[1][0]; dst_im[1][4] = mb_out_im[1][0];
    dst_re[1][5] = mb_out_re[1][1]; dst_im[1][5] = mb_out_im[1][1];
    dst_re[1][6] = mb_out_re[1][2]; dst_im[1][6] = mb_out_im[1][2];
    dst_re[1][7] = mb_out_re[1][3]; dst_im[1][7] = mb_out_im[1][3];
    dst_re[2][4] = mb_out_re[2][0]; dst_im[2][4] = mb_out_im[2][0];
    dst_re[2][5] = mb_out_re[2][1]; dst_im[2][5] = mb_out_im[2][1];
    dst_re[2][6] = mb_out_re[2][2]; dst_im[2][6] = mb_out_im[2][2];
    dst_re[2][7] = mb_out_re[2][3]; dst_im[2][7] = mb_out_im[2][3];
    dst_re[3][4] = mb_out_re[3][0]; dst_im[3][4] = mb_out_im[3][0];
    dst_re[3][5] = mb_out_re[3][1]; dst_im[3][5] = mb_out_im[3][1];
    dst_re[3][6] = mb_out_re[3][2]; dst_im[3][6] = mb_out_im[3][2];
    dst_re[3][7] = mb_out_re[3][3]; dst_im[3][7] = mb_out_im[3][3];

    dst_re[4][0] = mb_out_re[0][0]; dst_im[4][0] = VECT_NEG(mb_out_im[0][0]);
    dst_re[4][1] = mb_out_re[1][0]; dst_im[4][1] = VECT_NEG(mb_out_im[1][0]);
    dst_re[4][2] = mb_out_re[2][0]; dst_im[4][2] = VECT_NEG(mb_out_im[2][0]);
    dst_re[4][3] = mb_out_re[3][0]; dst_im[4][3] = VECT_NEG(mb_out_im[3][0]);
    dst_re[5][0] = mb_out_re[0][1]; dst_im[5][0] = VECT_NEG(mb_out_im[0][1]);
    dst_re[5][1] = mb_out_re[1][1]; dst_im[5][1] = VECT_NEG(mb_out_im[1][1]);
    dst_re[5][2] = mb_out_re[2][1]; dst_im[5][2] = VECT_NEG(mb_out_im[2][1]);
    dst_re[5][3] = mb_out_re[3][1]; dst_im[5][3] = VECT_NEG(mb_out_im[3][1]);
    dst_re[6][0] = mb_out_re[0][2]; dst_im[6][0] = VECT_NEG(mb_out_im[0][2]);
    dst_re[6][1] = mb_out_re[1][2]; dst_im[6][1] = VECT_NEG(mb_out_im[1][2]);
    dst_re[6][2] = mb_out_re[2][2]; dst_im[6][2] = VECT_NEG(mb_out_im[2][2]);
    dst_re[6][3] = mb_out_re[3][2]; dst_im[6][3] = VECT_NEG(mb_out_im[3][2]);
    dst_re[7][0] = mb_out_re[0][3]; dst_im[7][0] = VECT_NEG(mb_out_im[0][3]);
    dst_re[7][1] = mb_out_re[1][3]; dst_im[7][1] = VECT_NEG(mb_out_im[1][3]);
    dst_re[7][2] = mb_out_re[2][3]; dst_im[7][2] = VECT_NEG(mb_out_im[2][3]);
    dst_re[7][3] = mb_out_re[3][3]; dst_im[7][3] = VECT_NEG(mb_out_im[3][3]);

    dst_re[4][4] = md_out_re[0][0]; dst_im[4][4] = md_out_im[0][0];
    dst_re[4][5] = md_out_re[0][1]; dst_im[4][5] = md_out_im[0][1];
    dst_re[4][6] = md_out_re[0][2]; dst_im[4][6] = md_out_im[0][2];
    dst_re[4][7] = md_out_re[0][3]; dst_im[4][7] = md_out_im[0][3];
    dst_re[5][4] = md_out_re[1][0]; dst_im[5][4] = md_out_im[1][0];
    dst_re[5][5] = md_out_re[1][1]; dst_im[5][5] = md_out_im[1][1];
    dst_re[5][6] = md_out_re[1][2]; dst_im[5][6] = md_out_im[1][2];
    dst_re[5][7] = md_out_re[1][3]; dst_im[5][7] = md_out_im[1][3];
    dst_re[6][4] = md_out_re[2][0]; dst_im[6][4] = md_out_im[2][0];
    dst_re[6][5] = md_out_re[2][1]; dst_im[6][5] = md_out_im[2][1];
    dst_re[6][6] = md_out_re[2][2]; dst_im[6][6] = md_out_im[2][2];
    dst_re[6][7] = md_out_re[2][3]; dst_im[6][7] = md_out_im[2][3];
    dst_re[7][4] = md_out_re[3][0]; dst_im[7][4] = md_out_im[3][0];
    dst_re[7][5] = md_out_re[3][1]; dst_im[7][5] = md_out_im[3][1];
    dst_re[7][6] = md_out_re[3][2]; dst_im[7][6] = md_out_im[3][2];
    dst_re[7][7] = md_out_re[3][3]; dst_im[7][7] = md_out_im[3][3];

}

// Stores result in transposed form
inline void InverseMatrix8x8HTranspose(VECT_T src_re[8][8], VECT_T src_im[8][8],
                                VECT_T dst_re[8][8], VECT_T dst_im[8][8])
{
    VECT_T ma_re[4][4];
    VECT_T ma_im[4][4];
    VECT_T mb_re[4][4];
    VECT_T mb_im[4][4];
    VECT_T mc_re[4][4];
    VECT_T mc_im[4][4];
    VECT_T md_re[4][4];
    VECT_T md_im[4][4];

    VECT_T ma_out_re[4][4];
    VECT_T ma_out_im[4][4];
    VECT_T mb_out_re[4][4];
    VECT_T mb_out_im[4][4];
    VECT_T md_out_re[4][4];
    VECT_T md_out_im[4][4];

    VECT_T ma_inv_re[4][4];
    VECT_T ma_inv_im[4][4];

    VECT_T ma_temp_re[4][4];
    VECT_T ma_temp_im[4][4];

    VECT_T ma_temp2_re[4][4];
    VECT_T ma_temp2_im[4][4];

    VECT_T ma_cab_re[4][4];
    VECT_T ma_cab_im[4][4];

    VECT_T ma_ab_re[4][4];
    VECT_T ma_ab_im[4][4];

    VECT_T ma_ca_re[4][4];
    VECT_T ma_ca_im[4][4];

    // int ii;
    // int ij;

    // Create sub-matrices
    ma_re[0][0] = src_re[0][0]; ma_im[0][0] = src_im[0][0];
    ma_re[0][1] = src_re[0][1]; ma_im[0][1] = src_im[0][1];
    ma_re[0][2] = src_re[0][2]; ma_im[0][2] = src_im[0][2];
    ma_re[0][3] = src_re[0][3]; ma_im[0][3] = src_im[0][3];
    ma_re[1][0] = src_re[1][0]; ma_im[1][0] = src_im[1][0];
    ma_re[1][1] = src_re[1][1]; ma_im[1][1] = src_im[1][1];
    ma_re[1][2] = src_re[1][2]; ma_im[1][2] = src_im[1][2];
    ma_re[1][3] = src_re[1][3]; ma_im[1][3] = src_im[1][3];
    ma_re[2][0] = src_re[2][0]; ma_im[2][0] = src_im[2][0];
    ma_re[2][1] = src_re[2][1]; ma_im[2][1] = src_im[2][1];
    ma_re[2][2] = src_re[2][2]; ma_im[2][2] = src_im[2][2];
    ma_re[2][3] = src_re[2][3]; ma_im[2][3] = src_im[2][3];
    ma_re[3][0] = src_re[3][0]; ma_im[3][0] = src_im[3][0];
    ma_re[3][1] = src_re[3][1]; ma_im[3][1] = src_im[3][1];
    ma_re[3][2] = src_re[3][2]; ma_im[3][2] = src_im[3][2];
    ma_re[3][3] = src_re[3][3]; ma_im[3][3] = src_im[3][3];

    mb_re[0][0] = src_re[0][4]; mb_im[0][0] = src_im[0][4];
    mb_re[0][1] = src_re[0][5]; mb_im[0][1] = src_im[0][5];
    mb_re[0][2] = src_re[0][6]; mb_im[0][2] = src_im[0][6];
    mb_re[0][3] = src_re[0][7]; mb_im[0][3] = src_im[0][7];
    mb_re[1][0] = src_re[1][4]; mb_im[1][0] = src_im[1][4];
    mb_re[1][1] = src_re[1][5]; mb_im[1][1] = src_im[1][5];
    mb_re[1][2] = src_re[1][6]; mb_im[1][2] = src_im[1][6];
    mb_re[1][3] = src_re[1][7]; mb_im[1][3] = src_im[1][7];
    mb_re[2][0] = src_re[2][4]; mb_im[2][0] = src_im[2][4];
    mb_re[2][1] = src_re[2][5]; mb_im[2][1] = src_im[2][5];
    mb_re[2][2] = src_re[2][6]; mb_im[2][2] = src_im[2][6];
    mb_re[2][3] = src_re[2][7]; mb_im[2][3] = src_im[2][7];
    mb_re[3][0] = src_re[3][4]; mb_im[3][0] = src_im[3][4];
    mb_re[3][1] = src_re[3][5]; mb_im[3][1] = src_im[3][5];
    mb_re[3][2] = src_re[3][6]; mb_im[3][2] = src_im[3][6];
    mb_re[3][3] = src_re[3][7]; mb_im[3][3] = src_im[3][7];

    mc_re[0][0] = src_re[4][0]; mc_im[0][0] = src_im[4][0];
    mc_re[0][1] = src_re[4][1]; mc_im[0][1] = src_im[4][1];
    mc_re[0][2] = src_re[4][2]; mc_im[0][2] = src_im[4][2];
    mc_re[0][3] = src_re[4][3]; mc_im[0][3] = src_im[4][3];
    mc_re[1][0] = src_re[5][0]; mc_im[1][0] = src_im[5][0];
    mc_re[1][1] = src_re[5][1]; mc_im[1][1] = src_im[5][1];
    mc_re[1][2] = src_re[5][2]; mc_im[1][2] = src_im[5][2];
    mc_re[1][3] = src_re[5][3]; mc_im[1][3] = src_im[5][3];
    mc_re[2][0] = src_re[6][0]; mc_im[2][0] = src_im[6][0];
    mc_re[2][1] = src_re[6][1]; mc_im[2][1] = src_im[6][1];
    mc_re[2][2] = src_re[6][2]; mc_im[2][2] = src_im[6][2];
    mc_re[2][3] = src_re[6][3]; mc_im[2][3] = src_im[6][3];
    mc_re[3][0] = src_re[7][0]; mc_im[3][0] = src_im[7][0];
    mc_re[3][1] = src_re[7][1]; mc_im[3][1] = src_im[7][1];
    mc_re[3][2] = src_re[7][2]; mc_im[3][2] = src_im[7][2];
    mc_re[3][3] = src_re[7][3]; mc_im[3][3] = src_im[7][3];

    md_re[0][0] = src_re[4][4]; md_im[0][0] = src_im[4][4];
    md_re[0][1] = src_re[4][5]; md_im[0][1] = src_im[4][5];
    md_re[0][2] = src_re[4][6]; md_im[0][2] = src_im[4][6];
    md_re[0][3] = src_re[4][7]; md_im[0][3] = src_im[4][7];
    md_re[1][0] = src_re[5][4]; md_im[1][0] = src_im[5][4];
    md_re[1][1] = src_re[5][5]; md_im[1][1] = src_im[5][5];
    md_re[1][2] = src_re[5][6]; md_im[1][2] = src_im[5][6];
    md_re[1][3] = src_re[5][7]; md_im[1][3] = src_im[5][7];
    md_re[2][0] = src_re[6][4]; md_im[2][0] = src_im[6][4];
    md_re[2][1] = src_re[6][5]; md_im[2][1] = src_im[6][5];
    md_re[2][2] = src_re[6][6]; md_im[2][2] = src_im[6][6];
    md_re[2][3] = src_re[6][7]; md_im[2][3] = src_im[6][7];
    md_re[3][0] = src_re[7][4]; md_im[3][0] = src_im[7][4];
    md_re[3][1] = src_re[7][5]; md_im[3][1] = src_im[7][5];
    md_re[3][2] = src_re[7][6]; md_im[3][2] = src_im[7][6];
    md_re[3][3] = src_re[7][7]; md_im[3][3] = src_im[7][7];

    // Get inv(A)
    InverseMatrix4x4H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                4, 4, 4);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 4, 4, 4);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              4, 4);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix4x4H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                4, 4, 4);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 4, 4, 4);
    // [B]
    MATRIXSUB(g_mzero4x4, g_mzero4x4, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 4, 4);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 4, 4, 4);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 4, 4);

    // Put results back
    dst_re[0][0] = ma_out_re[0][0]; dst_im[0][0] = ma_out_im[0][0];
    dst_re[1][0] = ma_out_re[0][1]; dst_im[1][0] = ma_out_im[0][1];
    dst_re[2][0] = ma_out_re[0][2]; dst_im[2][0] = ma_out_im[0][2];
    dst_re[3][0] = ma_out_re[0][3]; dst_im[3][0] = ma_out_im[0][3];
    dst_re[0][1] = ma_out_re[1][0]; dst_im[0][1] = ma_out_im[1][0];
    dst_re[1][1] = ma_out_re[1][1]; dst_im[1][1] = ma_out_im[1][1];
    dst_re[2][1] = ma_out_re[1][2]; dst_im[2][1] = ma_out_im[1][2];
    dst_re[3][1] = ma_out_re[1][3]; dst_im[3][1] = ma_out_im[1][3];
    dst_re[0][2] = ma_out_re[2][0]; dst_im[0][2] = ma_out_im[2][0];
    dst_re[1][2] = ma_out_re[2][1]; dst_im[1][2] = ma_out_im[2][1];
    dst_re[2][2] = ma_out_re[2][2]; dst_im[2][2] = ma_out_im[2][2];
    dst_re[3][2] = ma_out_re[2][3]; dst_im[3][2] = ma_out_im[2][3];
    dst_re[0][3] = ma_out_re[3][0]; dst_im[0][3] = ma_out_im[3][0];
    dst_re[1][3] = ma_out_re[3][1]; dst_im[1][3] = ma_out_im[3][1];
    dst_re[2][3] = ma_out_re[3][2]; dst_im[2][3] = ma_out_im[3][2];
    dst_re[3][3] = ma_out_re[3][3]; dst_im[3][3] = ma_out_im[3][3];

    dst_re[4][0] = mb_out_re[0][0]; dst_im[4][0] = mb_out_im[0][0];
    dst_re[5][0] = mb_out_re[0][1]; dst_im[5][0] = mb_out_im[0][1];
    dst_re[6][0] = mb_out_re[0][2]; dst_im[6][0] = mb_out_im[0][2];
    dst_re[7][0] = mb_out_re[0][3]; dst_im[7][0] = mb_out_im[0][3];
    dst_re[4][1] = mb_out_re[1][0]; dst_im[4][1] = mb_out_im[1][0];
    dst_re[5][1] = mb_out_re[1][1]; dst_im[5][1] = mb_out_im[1][1];
    dst_re[6][1] = mb_out_re[1][2]; dst_im[6][1] = mb_out_im[1][2];
    dst_re[7][1] = mb_out_re[1][3]; dst_im[7][1] = mb_out_im[1][3];
    dst_re[4][2] = mb_out_re[2][0]; dst_im[4][2] = mb_out_im[2][0];
    dst_re[5][2] = mb_out_re[2][1]; dst_im[5][2] = mb_out_im[2][1];
    dst_re[6][2] = mb_out_re[2][2]; dst_im[6][2] = mb_out_im[2][2];
    dst_re[7][2] = mb_out_re[2][3]; dst_im[7][2] = mb_out_im[2][3];
    dst_re[4][3] = mb_out_re[3][0]; dst_im[4][3] = mb_out_im[3][0];
    dst_re[5][3] = mb_out_re[3][1]; dst_im[5][3] = mb_out_im[3][1];
    dst_re[6][3] = mb_out_re[3][2]; dst_im[6][3] = mb_out_im[3][2];
    dst_re[7][3] = mb_out_re[3][3]; dst_im[7][3] = mb_out_im[3][3];

    dst_re[0][4] = mb_out_re[0][0]; dst_im[0][4] = VECT_NEG(mb_out_im[0][0]);
    dst_re[1][4] = mb_out_re[1][0]; dst_im[1][4] = VECT_NEG(mb_out_im[1][0]);
    dst_re[2][4] = mb_out_re[2][0]; dst_im[2][4] = VECT_NEG(mb_out_im[2][0]);
    dst_re[3][4] = mb_out_re[3][0]; dst_im[3][4] = VECT_NEG(mb_out_im[3][0]);
    dst_re[0][5] = mb_out_re[0][1]; dst_im[0][5] = VECT_NEG(mb_out_im[0][1]);
    dst_re[1][5] = mb_out_re[1][1]; dst_im[1][5] = VECT_NEG(mb_out_im[1][1]);
    dst_re[2][5] = mb_out_re[2][1]; dst_im[2][5] = VECT_NEG(mb_out_im[2][1]);
    dst_re[3][5] = mb_out_re[3][1]; dst_im[3][5] = VECT_NEG(mb_out_im[3][1]);
    dst_re[0][6] = mb_out_re[0][2]; dst_im[0][6] = VECT_NEG(mb_out_im[0][2]);
    dst_re[1][6] = mb_out_re[1][2]; dst_im[1][6] = VECT_NEG(mb_out_im[1][2]);
    dst_re[2][6] = mb_out_re[2][2]; dst_im[2][6] = VECT_NEG(mb_out_im[2][2]);
    dst_re[3][6] = mb_out_re[3][2]; dst_im[3][6] = VECT_NEG(mb_out_im[3][2]);
    dst_re[0][7] = mb_out_re[0][3]; dst_im[0][7] = VECT_NEG(mb_out_im[0][3]);
    dst_re[1][7] = mb_out_re[1][3]; dst_im[1][7] = VECT_NEG(mb_out_im[1][3]);
    dst_re[2][7] = mb_out_re[2][3]; dst_im[2][7] = VECT_NEG(mb_out_im[2][3]);
    dst_re[3][7] = mb_out_re[3][3]; dst_im[3][7] = VECT_NEG(mb_out_im[3][3]);

    dst_re[4][4] = md_out_re[0][0]; dst_im[4][4] = md_out_im[0][0];
    dst_re[5][4] = md_out_re[0][1]; dst_im[5][4] = md_out_im[0][1];
    dst_re[6][4] = md_out_re[0][2]; dst_im[6][4] = md_out_im[0][2];
    dst_re[7][4] = md_out_re[0][3]; dst_im[7][4] = md_out_im[0][3];
    dst_re[4][5] = md_out_re[1][0]; dst_im[4][5] = md_out_im[1][0];
    dst_re[5][5] = md_out_re[1][1]; dst_im[5][5] = md_out_im[1][1];
    dst_re[6][5] = md_out_re[1][2]; dst_im[6][5] = md_out_im[1][2];
    dst_re[7][5] = md_out_re[1][3]; dst_im[7][5] = md_out_im[1][3];
    dst_re[4][6] = md_out_re[2][0]; dst_im[4][6] = md_out_im[2][0];
    dst_re[5][6] = md_out_re[2][1]; dst_im[5][6] = md_out_im[2][1];
    dst_re[6][6] = md_out_re[2][2]; dst_im[6][6] = md_out_im[2][2];
    dst_re[7][6] = md_out_re[2][3]; dst_im[7][6] = md_out_im[2][3];
    dst_re[4][7] = md_out_re[3][0]; dst_im[4][7] = md_out_im[3][0];
    dst_re[5][7] = md_out_re[3][1]; dst_im[5][7] = md_out_im[3][1];
    dst_re[6][7] = md_out_re[3][2]; dst_im[6][7] = md_out_im[3][2];
    dst_re[7][7] = md_out_re[3][3]; dst_im[7][7] = md_out_im[3][3];

}

/*void InverseMatrix16x16(VECT_T src_re[16][16], VECT_T src_im[16][16],
                      VECT_T dst_re[16][16], VECT_T dst_im[16][16])
{
    VECT_T ma_re[8][8];
    VECT_T ma_im[8][8];
    VECT_T mb_re[8][8];
    VECT_T mb_im[8][8];
    VECT_T mc_re[8][8];
    VECT_T mc_im[8][8];
    VECT_T md_re[8][8];
    VECT_T md_im[8][8];

    VECT_T ma_out_re[8][8];
    VECT_T ma_out_im[8][8];
    VECT_T mb_out_re[8][8];
    VECT_T mb_out_im[8][8];
    VECT_T mc_out_re[8][8];
    VECT_T mc_out_im[8][8];
    VECT_T md_out_re[8][8];
    VECT_T md_out_im[8][8];

    VECT_T ma_inv_re[8][8];
    VECT_T ma_inv_im[8][8];

    VECT_T ma_temp_re[8][8];
    VECT_T ma_temp_im[8][8];

    VECT_T ma_temp2_re[8][8];
    VECT_T ma_temp2_im[8][8];

    VECT_T ma_cab_re[8][8];
    VECT_T ma_cab_im[8][8];

    VECT_T ma_ab_re[8][8];
    VECT_T ma_ab_im[8][8];

    VECT_T ma_ca_re[8][8];
    VECT_T ma_ca_im[8][8];

    int ii;
    int ij;

    // Create sub-matrices
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            ma_re[ii][ij] = src_re[ii][ij];
            ma_im[ii][ij] = src_im[ii][ij];

            mb_re[ii][ij] = src_re[ii][ij + 8];
            mb_im[ii][ij] = src_im[ii][ij + 8];

            mc_re[ii][ij] = src_re[ii + 8][ij];
            mc_im[ii][ij] = src_im[ii + 8][ij];

            md_re[ii][ij] = src_re[ii + 8][ij + 8];
            md_im[ii][ij] = src_im[ii + 8][ij + 8];
        }
    }

    // Get inv(A)
    InverseMatrix8x8(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMUL(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
              8, 8, 8);

    // C*inv(A)*B
    MATRIXMUL(ma_ca_re, ma_ca_im, mb_re, mb_im, ma_temp2_re, ma_temp2_im,
              8, 8, 8);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              8, 8);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix8x8(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMUL(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
              8, 8, 8);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMUL(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
              ma_temp_re, ma_temp_im, 8, 8, 8);
    // [B]
    MATRIXSUB(g_mzero8x8, g_mzero8x8, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 8, 8);

    // (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(md_out_re, md_out_im, ma_ca_re, ma_ca_im,
              ma_temp_re, ma_temp_im, 8, 8, 8);

    // [C]
    MATRIXSUB(g_mzero8x8, g_mzero8x8, ma_temp_re, ma_temp_im,
              mc_out_re, mc_out_im, 8, 8);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMUL(ma_ab_re, ma_ab_im, ma_temp_re, ma_temp_im,
              ma_temp2_re, ma_temp2_im, 8, 8, 8);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 8, 8);

    // Put results back
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ii][ij] = ma_out_re[ii][ij];
            dst_im[ii][ij] = ma_out_im[ii][ij];

            dst_re[ii][ij + 8] = mb_out_re[ii][ij];
            dst_im[ii][ij + 8] = mb_out_im[ii][ij];

            dst_re[ii + 8][ij] = mc_out_re[ii][ij];
            dst_im[ii + 8][ij] = mc_out_im[ii][ij];

            dst_re[ii + 8][ij + 8] = md_out_re[ii][ij];
            dst_im[ii + 8][ij + 8] = md_out_im[ii][ij];
        }
    }
}
*/

/*void InverseMatrix16x16H(VECT_T src_re[16][16], VECT_T src_im[16][16],
                       VECT_T dst_re[16][16], VECT_T dst_im[16][16])
{
    VECT_T ma_re[8][8];
    VECT_T ma_im[8][8];
    VECT_T mb_re[8][8];
    VECT_T mb_im[8][8];
    VECT_T mc_re[8][8];
    VECT_T mc_im[8][8];
    VECT_T md_re[8][8];
    VECT_T md_im[8][8];

    VECT_T ma_out_re[8][8];
    VECT_T ma_out_im[8][8];
    VECT_T mb_out_re[8][8];
    VECT_T mb_out_im[8][8];
    VECT_T md_out_re[8][8];
    VECT_T md_out_im[8][8];

    VECT_T ma_inv_re[8][8];
    VECT_T ma_inv_im[8][8];

    VECT_T ma_temp_re[8][8];
    VECT_T ma_temp_im[8][8];

    VECT_T ma_temp2_re[8][8];
    VECT_T ma_temp2_im[8][8];

    VECT_T ma_cab_re[8][8];
    VECT_T ma_cab_im[8][8];

    VECT_T ma_ab_re[8][8];
    VECT_T ma_ab_im[8][8];

    VECT_T ma_ca_re[8][8];
    VECT_T ma_ca_im[8][8];

    int ii;
    int ij;

    // Create sub-matrices
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            ma_re[ii][ij] = src_re[ii][ij];
            ma_im[ii][ij] = src_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            mb_re[ii][ij] = src_re[ii][ij + 8];
            mb_im[ii][ij] = src_im[ii][ij + 8];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            mc_re[ii][ij] = src_re[ii + 8][ij];
            mc_im[ii][ij] = src_im[ii + 8][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            md_re[ii][ij] = src_re[ii + 8][ij + 8];
            md_im[ii][ij] = src_im[ii + 8][ij + 8];
        }
    }

    // Get inv(A)
    InverseMatrix8x8H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                8, 8, 8);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 8, 8, 8);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              8, 8);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix8x8H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                8, 8, 8);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 8, 8, 8);
    // [B]
    MATRIXSUB(g_mzero8x8, g_mzero8x8, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 8, 8);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 8, 8, 8);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 8, 8);

    // Put results back
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ii][ij] = ma_out_re[ii][ij];
            dst_im[ii][ij] = ma_out_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ii][ij + 8] = mb_out_re[ii][ij];
            dst_im[ii][ij + 8] = mb_out_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ii + 8][ij] = mb_out_re[ij][ii];
            dst_im[ii + 8][ij] = VECT_NEG(mb_out_im[ij][ii]);
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ii + 8][ij + 8] = md_out_re[ii][ij];
            dst_im[ii + 8][ij + 8] = md_out_im[ii][ij];
        }
    }
}
*/

// Stores result in transposed form
inline void InverseMatrix16x16HTranspose(VECT_T src_re[16][16], VECT_T src_im[16][16],
								  VECT_T dst_re[16][16], VECT_T dst_im[16][16])
{
    VECT_T ma_re[8][8];
    VECT_T ma_im[8][8];
    VECT_T mb_re[8][8];
    VECT_T mb_im[8][8];
    VECT_T mc_re[8][8];
    VECT_T mc_im[8][8];
    VECT_T md_re[8][8];
    VECT_T md_im[8][8];

    VECT_T ma_out_re[8][8];
    VECT_T ma_out_im[8][8];
    VECT_T mb_out_re[8][8];
    VECT_T mb_out_im[8][8];
    VECT_T md_out_re[8][8];
    VECT_T md_out_im[8][8];

    VECT_T ma_inv_re[8][8];
    VECT_T ma_inv_im[8][8];

    VECT_T ma_temp_re[8][8];
    VECT_T ma_temp_im[8][8];

    VECT_T ma_temp2_re[8][8];
    VECT_T ma_temp2_im[8][8];

    VECT_T ma_cab_re[8][8];
    VECT_T ma_cab_im[8][8];

    VECT_T ma_ab_re[8][8];
    VECT_T ma_ab_im[8][8];

    VECT_T ma_ca_re[8][8];
    VECT_T ma_ca_im[8][8];

    int ii;
    int ij;

    // Create sub-matrices
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            ma_re[ii][ij] = src_re[ii][ij];
            ma_im[ii][ij] = src_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            mb_re[ii][ij] = src_re[ii][ij + 8];
            mb_im[ii][ij] = src_im[ii][ij + 8];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            mc_re[ii][ij] = src_re[ii + 8][ij];
            mc_im[ii][ij] = src_im[ii + 8][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            md_re[ii][ij] = src_re[ii + 8][ij + 8];
            md_im[ii][ij] = src_im[ii + 8][ij + 8];
        }
    }

    // Get inv(A)
    InverseMatrix8x8H(ma_re, ma_im, ma_inv_re, ma_inv_im);

    // C*inv(A)
    MATRIXMULNH(mc_re, mc_im, ma_inv_re, ma_inv_im, ma_ca_re, ma_ca_im,
                8, 8, 8);

    // C*inv(A)*B
    MATRIXMULTHT(ma_ca_re, ma_ca_im, mb_re, mb_im,
                 ma_temp2_re, ma_temp2_im, 8, 8, 8);

    // D - C*inv(A)*B
    MATRIXSUB(md_re, md_im, ma_temp2_re, ma_temp2_im, ma_cab_re, ma_cab_im,
              8, 8);

    // inv(D - C*inv(A)*B) -> [D]
    InverseMatrix8x8H(ma_cab_re, ma_cab_im, md_out_re, md_out_im);

    // inv(A)*B
    MATRIXMULHN(ma_inv_re, ma_inv_im, mb_re, mb_im, ma_ab_re, ma_ab_im,
                8, 8, 8);

    // (inv(A)*B) * (inv(D - C*inv(A)*B))
    MATRIXMULNH(ma_ab_re, ma_ab_im, md_out_re, md_out_im,
                ma_temp_re, ma_temp_im, 8, 8, 8);
    // [B]
    MATRIXSUB(g_mzero8x8, g_mzero8x8, ma_temp_re, ma_temp_im,
              mb_out_re, mb_out_im, 8, 8);

    // inv(A)*B * (inv(D - C*inv(A)*B)) * C*inv(A)
    MATRIXMULTHT(ma_temp_re, ma_temp_im, ma_ca_re, ma_ca_im,
                 ma_temp2_re, ma_temp2_im, 8, 8, 8);

    // [A]
    MATRIXADD(ma_inv_re, ma_inv_im, ma_temp2_re, ma_temp2_im,
              ma_out_re, ma_out_im, 8, 8);

    // Put results back
    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ij][ii] = ma_out_re[ii][ij];
            dst_im[ij][ii] = ma_out_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ij + 8][ii] = mb_out_re[ii][ij];
            dst_im[ij + 8][ii] = mb_out_im[ii][ij];
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ij][ii + 8] = mb_out_re[ij][ii];
            dst_im[ij][ii + 8] = VECT_NEG(mb_out_im[ij][ii]);
        }
    }

    for (ii = 0; ii < 8; ii++) {
        for (ij = 0; ij < 8; ij++) {
            dst_re[ij + 8][ii + 8] = md_out_re[ii][ij];
            dst_im[ij + 8][ii + 8] = md_out_im[ii][ij];
        }
    }
}

/*void InverseMatrix2x2FP(VECT_I src[2][2], VECT_T dst_re[2][2],
                        VECT_T dst_im[2][2])
{
    VECT_T src_re[2][2];
    VECT_T src_im[2][2];

    MATRIXCVT16FXTO32FP(src, src_re, src_im, 2, 2, gFactor1);

    InverseMatrix2x2(src_re, src_im, dst_re, dst_im);
}

void InverseMatrix2x2HFP(VECT_I src[2][2], VECT_T dst_re[2][2],
                         VECT_T dst_im[2][2])
{
    VECT_T src_re[2][2];
    VECT_T src_im[2][2];

    MATRIXCVT16FXTO32FPH(src, src_re, src_im, 2, gFactor1);

    InverseMatrix2x2H(src_re, src_im, dst_re, dst_im);
}

void InverseMatrix4x4FP(VECT_I src[4][4], VECT_T dst_re[4][4],
                        VECT_T dst_im[4][4])
{
    VECT_T src_re[4][4];
    VECT_T src_im[4][4];

    MATRIXCVT16FXTO32FP(src, src_re, src_im, 4, 4, gFactor1);

    InverseMatrix4x4(src_re, src_im, dst_re, dst_im);
}

void InverseMatrix4x4HFP(VECT_I src[4][4], VECT_T dst_re[4][4],
                         VECT_T dst_im[4][4])
{
    VECT_T src_re[4][4];
    VECT_T src_im[4][4];

    MATRIXCVT16FXTO32FPH(src, src_re, src_im, 4, gFactor1);

    InverseMatrix4x4H(src_re, src_im, dst_re, dst_im);
}

void InverseMatrix8x8FP(VECT_I src[8][8], VECT_T dst_re[8][8],
                        VECT_T dst_im[8][8])
{
    VECT_T src_re[8][8];
    VECT_T src_im[8][8];

    MATRIXCVT16FXTO32FP(src, src_re, src_im, 8, 8, gFactor1);

    InverseMatrix8x8(src_re, src_im, dst_re, dst_im);
}

void InverseMatrix8x8HFP(VECT_I src[8][8], VECT_T dst_re[8][8],
                         VECT_T dst_im[8][8])
{
    VECT_T src_re[8][8];
    VECT_T src_im[8][8];

    MATRIXCVT16FXTO32FPH(src, src_re, src_im, 8, gFactor1);

    InverseMatrix8x8H(src_re, src_im, dst_re, dst_im);
}
*/

#endif /* _MATRIX_INVERT_INC_ */
