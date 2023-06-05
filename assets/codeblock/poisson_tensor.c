void tabulate_tensor_integral_c2833a6d26cdde02c09f8cbb4b3bf92016ec9e7e(double* restrict A,
                                    const double* restrict w,
                                    const double* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{
  // Quadrature rules
  static const double weights_083[1] = { 0.5 };
  // Precomputed values of basis functions and precomputations
  // FE* dimensions: [permutation][entities][points][dofs]
  static const double FE3_C0_D01_Q083[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
  static const double FE4_C0_D10_Q083[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
  // Quadrature loop independent computations for quadrature rule 083
  double J_c0 = 0.0;
  double J_c3 = 0.0;
  double J_c1 = 0.0;
  double J_c2 = 0.0;
  for (int ic = 0; ic < 3; ++ic)
  {
    J_c0 += coordinate_dofs[ic * 3] * FE4_C0_D10_Q083[0][0][0][ic];
    J_c3 += coordinate_dofs[ic * 3 + 1] * FE3_C0_D01_Q083[0][0][0][ic];
    J_c1 += coordinate_dofs[ic * 3] * FE3_C0_D01_Q083[0][0][0][ic];
    J_c2 += coordinate_dofs[ic * 3 + 1] * FE4_C0_D10_Q083[0][0][0][ic];
  }
  double sp_083[20];
  sp_083[0] = J_c0 * J_c3;
  sp_083[1] = J_c1 * J_c2;
  sp_083[2] = sp_083[0] + -1 * sp_083[1];
  sp_083[3] = J_c0 / sp_083[2];
  sp_083[4] = (-1 * J_c1) / sp_083[2];
  sp_083[5] = sp_083[3] * sp_083[3];
  sp_083[6] = sp_083[3] * sp_083[4];
  sp_083[7] = sp_083[4] * sp_083[4];
  sp_083[8] = J_c3 / sp_083[2];
  sp_083[9] = (-1 * J_c2) / sp_083[2];
  sp_083[10] = sp_083[9] * sp_083[9];
  sp_083[11] = sp_083[8] * sp_083[9];
  sp_083[12] = sp_083[8] * sp_083[8];
  sp_083[13] = sp_083[5] + sp_083[10];
  sp_083[14] = sp_083[6] + sp_083[11];
  sp_083[15] = sp_083[12] + sp_083[7];
  sp_083[16] = fabs(sp_083[2]);
  sp_083[17] = sp_083[13] * sp_083[16];
  sp_083[18] = sp_083[14] * sp_083[16];
  sp_083[19] = sp_083[15] * sp_083[16];
  for (int iq = 0; iq < 1; ++iq)
  {
    const double fw0 = sp_083[19] * weights_083[iq];
    const double fw1 = sp_083[18] * weights_083[iq];
    const double fw2 = sp_083[17] * weights_083[iq];
    double t0[3];
    double t1[3];
    for (int i = 0; i < 3; ++i)
    {
      t0[i] = fw0 * FE4_C0_D10_Q083[0][0][0][i] + fw1 * FE3_C0_D01_Q083[0][0][0][i];
      t1[i] = fw1 * FE4_C0_D10_Q083[0][0][0][i] + fw2 * FE3_C0_D01_Q083[0][0][0][i];
    }
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        A[3 * i + j] += FE4_C0_D10_Q083[0][0][0][j] * t0[i] + FE3_C0_D01_Q083[0][0][0][j] * t1[i];
  }
}