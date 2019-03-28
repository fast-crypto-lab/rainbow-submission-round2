# rainbow

Rainbow Signature system.


#usage:

0. (Optional) Edit rainbow_config.h for setting (V1,O1,O2).
   The default: (gf(16),32,32,32).

1. Make for 3 executables: rainbow-genkey, rainbow-sign, and rainbow-verify .

2. (generate key pairs)
  rainbow-genkey  pk_file_name  sk_file_name [selfdefined_random_file]

3. (sign a file)
  rainbow-sign  sk_file_name file_name_to_be_signed

   (Or redirect the signature to an output file)
  rainbow-sign  sk_file_name file_name_to_be_signed | tee signature_file_name


4. (verify a signature)
  rainbow-verify  pk_file_name  signature_file_name  file_name_to_be_signed



# (Optional) Use self-defined randomness:

1. Prepare a binary file and add the file name to the command line while executing.
   The contents of the file are the ``seed of randomness'' used.

