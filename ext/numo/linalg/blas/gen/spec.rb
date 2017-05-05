if /[cz]/ =~ blas_char
  def_id "real"
  def_id "imag"
end

case blas_char
when "c"
  real_char = "s"
when "z"
  real_char = "d"
end

if /[sd]/ =~ blas_char
  def_blas "?dot"
  def_blas "?nrm2"
  def_blas "?asum", "nrm2"
end
if "s" == blas_char
  def_blas "dsdot", "dot"
  def_blas "sdsdot"
end
if /[cz]/ =~ blas_char
  def_blas "?dotc", "dot", "?dotc_sub"
  def_blas "?dotu", "dot", "?dotu_sub"
  def_blas real_char+"?nrm2", "nrm2"
  def_blas real_char+"?asum", "nrm2"
end

def_blas "?swap"
def_blas "?copy"
def_blas "?axpy"

if /[sd]/ =~ blas_char
  def_blas "?rot"
  def_blas "?rotm"
else
  #def_blas blas_char+real_char+"rot", "rot"
end

def_blas "?scal"
if /[cz]/ =~ blas_char
  def_blas blas_char+real_char+"scal", "scal"
end

# level-2 blas
def_blas "?gemv"
def_blas "?trmv"
if /[cz]/ =~ blas_char
  def_blas "?hemv"
end
case blas_char
when /[sd]/
  def_blas "?ger"
  def_blas "?syr"
  def_blas "?syr2"
when /[cz]/
  def_blas "?gerc", "ger"
  def_blas "?geru", "ger"
  def_blas "?her",  "syr"
  def_blas "?her2", "syr2"
end

# level-3 blas
def_blas "?gemm"
def_blas "?symm"
if /[cz]/ =~ blas_char
  def_blas "?hemm", "symm"
end
def_blas "?syrk"
def_blas "?syr2k"
