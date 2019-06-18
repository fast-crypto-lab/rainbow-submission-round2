
rm_kats()
{
  echo $1
  echo "rm KATs..."
  rm $1/*.req $1/*.rsp
  return 0
}


### Main script start ###

echo " ================ rm KATs ================"

rm_kats Ia_Classic
rm_kats Ia_Cyclic
rm_kats Ia_CompressedCyclic
rm_kats IIIc_Classic
rm_kats IIIc_Cyclic
rm_kats IIIc_CompressedCyclic
rm_kats Vc_Classic
rm_kats Vc_Cyclic
rm_kats Vc_CompressedCyclic
