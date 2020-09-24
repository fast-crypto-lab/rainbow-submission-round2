# rainbow

Rainbow Signature system.


#usage:

0. (Optional) Edit the makefile to choose a specific variant to build.
   For example of building the `Ia_Classic' variant, 
   set the variable $PROJ_DIR = Ia_Classic in the makefile.
   It will build the codes in the `Ia_Classic' folder only.

   The default variant is `Ia_Classic' and all possiable variants are 
   `Ia_Classic',   `Ia_Circumzenithal',       `Ia_Compressed', 
   `IIIc_Classic', `IIIc_Circumzenithal',     `IIIc_Compressed',
   `Vc_Classic',   `Vc_Circumzenithal',   and `Vc_Compressed'.

   The $PROJ_DIR can be setted either in the makefile explictly or in command-line. 
   No code changes have to be made.


1. Make for 3 executables: rainbow-genkey, rainbow-sign, and rainbow-verify .

   Example 1: build the Ia_Compressed variant, type:
   make PROJ_DIR=Ia_Compressed

   Example 2: build the Vc_Classic variant, type:
   make PROJ_DIR=Vc_Classic


2. (generate key pairs)
  rainbow-genkey  pk_file_name  sk_file_name [selfdefined_random_file]

3. (sign a file)
  rainbow-sign  sk_file_name file_name_to_be_signed

   (Or redirect the signature to an output file)
  rainbow-sign  sk_file_name file_name_to_be_signed | tee signature_file_name


4. (verify a signature)
  rainbow-verify  pk_file_name  signature_file_name  file_name_to_be_signed


5. (generate KATs) :
  genKATs.sh
  execute `genKATs.sh'. The KATs will be copied to corresponding source folders.
