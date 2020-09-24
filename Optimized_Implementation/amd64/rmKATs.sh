
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
rm_kats Ia_Circumzenithal
rm_kats Ia_Compressed
rm_kats IIIc_Classic
rm_kats IIIc_Circumzenithal
rm_kats IIIc_Compressed
rm_kats Vc_Classic
rm_kats Vc_Circumzenithal
rm_kats Vc_Compressed
