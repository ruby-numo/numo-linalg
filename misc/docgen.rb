class Par
  def initialize(name)
    @name = name
    @docs = {}
    @docs["summary"] = @doc = []
  end

  def add(line)
    if @verb
      @doc << line
    end
  end

  def command(cmd,param)
    case cmd
    when "verbatim"
      @verb = true
    when "endverbatim"
      @verb = false
    when /^param(.*)$/
      inout = $1
      h = (@docs["param"] ||= {})
      h[param.upcase] = @doc = [inout]
    end
  end

  def finish
    @doc = nil
  end

  def [](key)
    @docs[key]
  end
end

class FuncDoc
  def initialize(fname)
    @name = File.basename(fname,'.f')
    @par = {}
  end
  attr_reader :name

  def add(line)
    case line
    when /^\\(\S+)(\s+(.*))?$/
      command($1,$3)
    else
      @buf.add(line) if @buf
    end
  end

  def command(cmd,param)
    case cmd
    when "par"
      @buf.finish if @buf
      @par[param] = @buf = Par.new(param)
    else
      @buf.command(cmd,param) if @buf
    end
  end

  def finish
    @buf.finish if @buf
    @buf = nil
  end

  def [](key)
    @par[key]
  end
end


def docread(dir)
  blas_dir = File.join(ENV['HOME'],dir)
  docs = {}
  Dir.glob(blas_dir+"/*.f").each do |src|
    dg = FuncDoc.new(src)
    open(src,'rt') do |f|
      f.each_line do |line|
        case line
        when /^\*> (.*)$/
          dg.add($1)
        end
      end
    end
    dg.finish
    if d = dg['Purpose:']
      #docs[dg.name] = d['summary']
      docs[dg.name] = h = {}
      h['summary'] = d['summary']
      h['param'] = d['param']
    end
  end
  docs
end

require 'pp'

def gendoc(outfile,dir)
  doc = docread(dir)
  open(outfile,"w"){|f| PP.pp(doc,f) }
end

gendoc("../ext/numo/linalg/lapack/gen/desc.rb","src/lapack-3.7.0/SRC")
gendoc("../ext/numo/linalg/blas/gen/desc.rb","src/lapack-3.7.0/BLAS/SRC")
