How to generate KATs

1. moving the working directory to 'Alternative_Implementation/avx2' and execute 'genKATs.sh'
  cd Alternative_Implementation/avx2
  . ./genKATs.sh

2. Collect KATs: execute 'cpKATs.sh'. It will collect all KATs in the 'KAT' directory 
  . ./cpKATs.sh

3. Move the KAT fiels in the 'KAT' directory to your target position.
  mv KAT ~/my_KAT

4. Clean the generated KAT files.
  . ./rmKATs.sh
