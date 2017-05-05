require_relative './erbpp2'

module BlasMethod
  def def_blas(meth, tmpl=nil, cblasfn=nil, **h)
    c = get(:blas_char)
    tmpl = meth.dup.gsub(/\?/,"") unless tmpl
    meth = meth.gsub(/\?/,c)
    h[:cblas_func] ||= (cblasfn && cblasfn.gsub(/\?/,c)) || meth

    case meth
    when /^s/
      result_dtype = "float"
      result_class = "numo_cSFloat"
    when /^d/
      result_dtype = "double"
      result_class = "numo_cDFloat"
    when /^c/
      result_dtype = "scomplex"
      result_class = "numo_cSComplex"
    when /^z/
      result_dtype = "dcomplex"
      result_class = "numo_cDComplex"
    end
    h[:result_dtype] ||= result_dtype
    h[:result_class] ||= result_class

    def_module_function(meth, tmpl, **h)
  end
end
