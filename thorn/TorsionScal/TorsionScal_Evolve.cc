void TorsionScal_Evolve(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  for (int k=1; k<cctk_lsh[2]-1; ++k)
  for (int j=1; j<cctk_lsh[1]-1; ++j)
  for (int i=1; i<cctk_lsh[0]-1; ++i)
  {
    const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    const double laplacian =
      (Torsion[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - 2.0*Torsion[index] + Torsion[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) / dx2 +
      (Torsion[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - 2.0*Torsion[index] + Torsion[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) / dy2 +
      (Torsion[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - 2.0*Torsion[index] + Torsion[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) / dz2;

    Torsion_rhs[index] = beta * laplacian
                       - lambda * (Torsion[index]*Torsion[index] - S_Pl*S_Pl) * Torsion[index];
  }
}
