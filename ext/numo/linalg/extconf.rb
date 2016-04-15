require 'rbconfig.rb'
require 'numo/narray'

#RbConfig::MAKEFILE_CONFIG["optflags"] = "-g3 -gdwarf-2"

require 'mkmf'

# configure options: --with-lapack-lib=path
dir_config("lapack")

# configure options: --with-blas-lib=path
dir_config("blas")

if with_config("openblas")
  # configure options: --with-openblas
  exit(-1) unless have_library("openblas")
elsif with_config("lapack")
  # configure options: --with-lapack
  exit(-1) unless have_library("blas") && have_library("lapack")
elsif have_library("openblas")
  # check openblas
elsif have_library("blas") && have_library("lapack")
  # check lapack
else
  puts "LAPACK library was not found."
  exit(-1)
end

#$CFLAGS="-g -O0"

$LOAD_PATH.each do |x|
  if File.exist? File.join(x,'numo/numo/narray.h')
    $INCFLAGS = "-I#{x}/numo " + $INCFLAGS
    break
  end
end

srcs = %w(
linalg
linalg_d
linalg_s
linalg_z
linalg_c
)
$objs = srcs.collect{|i| i+".o"}
fflags = ""

# # GNU FORTRAN v4
# if have_library("gfortran")
#   $defs.push "-fPIC -DGNU_FORTRAN"
#   fc = "gfortran"
#   if false # have_library("gomp")
#     fflags += "-fopenmp"
#   end
# #
# # GNU FORTRAN v3
# elsif have_library("g77")
#   $defs.push "-fPIC -DGNU_FORTRAN"
#   fc = "g77"
# elsif have_library('f2c')
#   $defs.push "-DF2C"
# else
#   puts "failed"
#   exit
# end

if !have_header('numo/narray.h')
  print <<EOL
  Header numo/narray.h was not found. Give pathname as follows:
  % ruby extconf.rb --with-narray-include=narray_h_dir
EOL
  exit(-1)
end

system "rm -f depend; erb depend.erb > depend"

create_makefile('numo/linalg')

# if $makefile_created
#   puts "appending extra install tasks to Makefile"
#   File.open("Makefile","a") do |w|
#     w.print <<EOL
# F77 = #{fc}
# FFLAGS = -O3 -fPIC -Iffte-5.0 #{fflags} -fomit-frame-pointer
# EOL
#   end
# end
