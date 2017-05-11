require_relative './erbpp2'

module LapackMethod
  def def_lpk(meth, tmpl=nil, fn=nil, **h)
    c = get(:blas_char)
    tmpl = meth.dup.gsub(/\?/,"") unless tmpl
    meth = meth.gsub(/\?/,c)
    h[:func_name] ||= (fn && fn.gsub(/\?/,c)) || meth

    case meth
    when /^s/
      result_dtype = "float"
      result_class = "numo_cSFloat"
      class_name = "SFloat"
      real_class_name = "SFloat"
      complex_class_name = "SComplex"
    when /^d/
      result_dtype = "double"
      result_class = "numo_cDFloat"
      class_name = "DFloat"
      real_class_name = "DFloat"
      complex_class_name = "DComplex"
    when /^c/
      result_dtype = "scomplex"
      result_class = "numo_cSComplex"
      class_name = "SComplex"
      real_class_name = "SFloat"
      complex_class_name = "SComplex"
    when /^z/
      result_dtype = "dcomplex"
      result_class = "numo_cDComplex"
      class_name = "DComplex"
      real_class_name = "DFloat"
      complex_class_name = "DComplex"
    end
    h[:result_dtype] ||= result_dtype
    h[:result_class] ||= result_class
    h[:class_name] ||= class_name
    h[:real_class_name] ||= real_class_name
    h[:complex_class_name] ||= complex_class_name

    case c
    when /c|z/
      h[:is_complex] ||= true
    else
      h[:is_complex] ||= false
    end

    def_module_function(meth, tmpl, **h)
  end
end
