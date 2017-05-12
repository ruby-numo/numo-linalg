def_id "order"
def_id "uplo"
def_id "jobu"
def_id "jobvt"
def_id "jobz"
def_id "jobvl"
def_id "jobvr"
def_id "trans"
def_id "rcond"
def_id "itype"

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

#--

def_lpk "?gesv"
def_lpk "?sysv", "gesv"
def_lpk "?posv", "gesv"
case blas_char
when /c|z/
  def_lpk "?hesv", "gesv"
end

def_lpk "?gesvd"
def_lpk "?gesdd", "gesvd"

def_lpk "?gels"
def_lpk "?gelss", "gels"
def_lpk "?gelsd", "gels"
def_lpk "?gelsy", "gels"

def_lpk "?geev"
def_lpk "?ggev"

case blas_char
when /c|z/
  def_lpk "?heev", "syev"
  def_lpk "?heevd", "syev"
  def_lpk "?hegv", "sygv"
  def_lpk "?hegvd", "sygv"
else
  def_lpk "?syev"
  def_lpk "?syevd", "syev"
  def_lpk "?sygv"
  def_lpk "?sygvd", "sygv"
end

def_lpk "?geqrf", "fact"
def_lpk "?gerqf", "fact"
def_lpk "?geqlf", "fact"
def_lpk "?gelqf", "fact"

def_lpk "?geqp3", "fact"
def_lpk "?tzrzf", "fact"

def_lpk "?getrf", "fact"
def_lpk "?sytrf", "fact"
def_lpk "?potrf", "fact"

def_lpk "?getri", "fact"
def_lpk "?sytri", "fact"
def_lpk "?potri", "fact"

case blas_char
when /c|z/
  def_lpk "?hetrf", "fact"
  def_lpk "?hetri", "fact"
end
