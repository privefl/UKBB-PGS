system("git clone https://github.com/getian107/PRScsx.git")

bigassertr::assert_dir("PRScsx/ldref")

library(runonce)
dl_ldref <- function(file, dir = "PRScsx/ldref") {
  tgz <- download_file(paste0("https://www.dropbox.com/s/", file, "?dl=1"),
                       dir = dir)
  untar(tgz, exdir = dir)
}
for (file in c("dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz",
               "fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz",
               "t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz",
               "nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz"))
  dl_ldref(file)

download_file("https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=1",
              dir = "PRScsx/ldref")
