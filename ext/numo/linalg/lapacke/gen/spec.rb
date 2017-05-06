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


def_lpk "?gesv"
def_lpk "?sysv", "gesv"
def_lpk "?posv", "gesv"
case blas_char
when /c|z/
  def_lpk "?hesv", "gesv"
end
