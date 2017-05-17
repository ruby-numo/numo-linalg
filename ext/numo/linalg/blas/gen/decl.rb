require_relative './erbpp2'

class DefLinalgFunction < DefModuleFunction
  def param(v,dim=nil,tp=nil,**h)
    tp ||= class_name
    inpl = h[:inplace] || ", inplace allowed"

    case dim
    when 2
      return "@param [#{tp}] #{v}  matrix (>=2-dimentional NArray#{inpl})."
    when 1
      return "@param [#{tp}] #{v}  vector (>=1-dimentional NArray#{inpl})."
    end

    optp = "String,Symbol"
    case v
    when /^order$/
      "@param [#{optp}] #{v}  if 'R': Row-major, if 'C': Column-major. (default='R')"
    when /^uplo$/
      "@param [#{optp}] #{v}  if 'U': Upper triangle, if 'L': Lower triangle. (default='U')"
    when /^side$/
      "@param [#{optp}] #{v}  if 'L': op(A)*B (left-side op), if 'R': B*op(A) (right-side op). (default='L')"
    when /^diag$/
      "@param [#{optp}] #{v}  if 'U': assumed to be unit triangular, if 'N': not assumed to be unit triangular. (default='U')"
    when /^trans(\w+)?$/
      a = dim || $1
      "@param [#{optp}] #{v}  if 'N': Not transpose #{a}, if 'T': Transpose #{a}. (default='N')"
    when /^job(\w+)?$/
      a = dim || $1
      "@param [#{optp}] #{v}  if 'V': Compute #{a}, if 'N': Not compute #{a} (default='V')"
    else
      raise "param not found: #{v}"
    end
  end
end

module Decl
  def decl(meth, tmpl=nil, fn=nil, **h)
    c = get(:blas_char)
    tmpl = meth.dup.gsub(/\?/,"") unless tmpl
    meth = meth.gsub(/\?/,c)
    h[:func_name] ||= (fn && fn.gsub(/\?/,c)) || meth
    h[:description] ||= get_desc(meth)

    case meth
    when /^s/
      result_dtype = "float"
      result_class = "numo_cSFloat"
      class_name = "Numo::SFloat"
      real_class_name = "Numo::SFloat"
      complex_class_name = "Numo::SComplex"
    when /^d/
      result_dtype = "double"
      result_class = "numo_cDFloat"
      class_name = "Numo::DFloat"
      real_class_name = "Numo::DFloat"
      complex_class_name = "Numo::DComplex"
    when /^c/
      result_dtype = "scomplex"
      result_class = "numo_cSComplex"
      class_name = "Numo::SComplex"
      real_class_name = "Numo::SFloat"
      complex_class_name = "Numo::SComplex"
    when /^z/
      result_dtype = "dcomplex"
      result_class = "numo_cDComplex"
      class_name = "Numo::DComplex"
      real_class_name = "Numo::DFloat"
      complex_class_name = "Numo::DComplex"
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

    DefLinalgFunction.new(self, tmpl, name:meth, **h)
  end

  def get_desc(meth)
    indent = false
    lines = []
    input = get(:module_desc)[meth]
    if input.nil?
      puts "no description: #{meth}"
      return ""
    end
    re_trim = /^(\s+)/=~input[0] ? /^#{$1}/ : nil
    input.each do |line|
      line.sub!(re_trim,"") if re_trim
      case line
      when /^\s*\d+\. /, /^\s*\(\d+\)/   # itemize
        lines << ""
        indent ||= ""
        line = indent+line
      when /^(\s+)/ #, /^(\s*)\w+\s*=/  # indent
        n_space = $1.size
        if !indent
          lines << ""
          indent = " " * [0,8-n_space].max
        end
        line = indent+line
      else  # standard
        if indent
          lines << ""
          indent = false
        end
      end
      if /^\s{8,}/ !~ line
        line.gsub!(/\*/,"\\*")
      end
      lines << line
    end
    lines.join("\n")+"\n"
  end
end
