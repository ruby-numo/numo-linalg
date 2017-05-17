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
def_id "norm"
def_id "axis"
def_id "keepdims"

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

# norm

decl "?lange"

# linear system

decl "?gesv"
decl "?sysv", "gesv"
decl "?posv", "gesv"
case blas_char
when /c|z/
  decl "?hesv", "gesv"
end

decl "?gesvd"
decl "?gesdd", "gesvd"

decl "?gels"
decl "?gelss", "gels"
decl "?gelsd", "gels"
decl "?gelsy", "gels"

# eigenvalue

decl "?geev"
decl "?ggev"

case blas_char
when /c|z/
  decl "?heev", "syev"
  decl "?heevd", "syev"
  decl "?hegv", "sygv"
  decl "?hegvd", "sygv"
else
  decl "?syev"
  decl "?syevd", "syev"
  decl "?sygv"
  decl "?sygvd", "sygv"
end

# factorize

decl "?geqrf", "fact"
decl "?gerqf", "fact"
decl "?geqlf", "fact"
decl "?gelqf", "fact"

decl "?geqp3", "fact"
decl "?tzrzf", "fact"

decl "?getrf", "fact"
decl "?sytrf", "fact"
decl "?potrf", "fact"

decl "?getri", "fact"
decl "?sytri", "fact"
decl "?potri", "fact"

case blas_char
when /c|z/
  decl "?hetrf", "fact"
  decl "?hetri", "fact"
end
