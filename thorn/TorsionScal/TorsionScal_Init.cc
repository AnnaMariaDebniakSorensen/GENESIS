void TorsionScal_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  for (int k=0; k<cctk_lsh[2]; ++k)
  for (int j=0; j<cctk_lsh[1]; ++j)
  for (int i=0; i<cctk_lsh[0]; ++i)
  {
    const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    Torsion[index] = S_Pl * exp(-r[index]);  // Initial radial profile
    Torsion_rhs[index] = 0.0;
  }
}
