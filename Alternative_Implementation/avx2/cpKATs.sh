
cp_kats()
{
  echo $1
  echo "cp KATs..."
  mkdir KAT/$1
  cp $1/*.req $1/*.rsp KAT/$1/
  return 0
}


### Main script start ###

echo " ================ cp KATs ================"
mkdir KAT

cp_kats Ia_Classic
cp_kats Ia_Circumzenithal
cp_kats Ia_Compressed
cp_kats IIIc_Classic
cp_kats IIIc_Circumzenithal
cp_kats IIIc_Compressed
cp_kats Vc_Classic
cp_kats Vc_Circumzenithal
cp_kats Vc_Compressed
