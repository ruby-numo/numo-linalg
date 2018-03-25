require 'mkmf'

def create_site_conf
  ldirs = [
    dir_config("mkl")[1],
    dir_config("openblas")[1],
    dir_config("atlas")[1],
    dir_config("blas")[1],
    dir_config("lapack")[1],
  ]
  bked = with_config("backend")

  require 'fiddle'

  message "creating lib/site_conf.rb\n"

  FileUtils.mkdir_p "lib"

  ext = detect_library_extension
  need_version = false
  if ext == 'so'
    begin
      Fiddle.dlopen "libm.so"
    rescue
      (5..7).each do |i|
        Fiddle.dlopen "libm.so.#{i}"
        need_version = true
        break
      rescue
      end
      if !need_version
        raise "failed to check whether dynamically linked shared object needs version suffix"
      end
    end
  end

  open("lib/site_conf.rb","w"){|f| f.write "
module Numo
  module Linalg

    BACKEND = #{bked.inspect}
    MKL_LIBPATH = #{ldirs[0].inspect}
    OPENBLAS_LIBPATH = #{ldirs[1].inspect}
    ATLAS_LIBPATH = #{ldirs[2].inspect}
    BLAS_LIBPATH = #{ldirs[3].inspect}
    LAPACK_LIBPATH = #{ldirs[4].inspect}

    module Loader
      EXT = '#{ext}'
      NEED_VERSION_SUFFIX = #{need_version}
    end

  end
end"
  }
end

def detect_library_extension
  case RbConfig::CONFIG['host_os']
  when /mswin|msys|mingw|cygwin/
    'dll'
  when /darwin|mac os/
    'dylib'
  else
    'so'
  end
end

def find_narray_h
  $LOAD_PATH.each do |x|
    if File.exist? File.join(x,'numo/numo/narray.h')
      $INCFLAGS = "-I#{x}/numo " + $INCFLAGS
      break
    end
  end
end

def find_libnarray_a
  $LOAD_PATH.each do |x|
    if File.exist? File.join(x,'numo/libnarray.a')
      $LDFLAGS = "-L#{x}/numo " + $LDFLAGS
      break
    end
  end
end

def create_depend
  require 'erb'
  message "creating depend\n"
  dep_path = File.join(Dir.pwd, "depend")
  File.open(dep_path, "w") do |dep|
    dep_erb_path = File.join(Dir.pwd, "depend.erb")
    File.open(dep_erb_path, "r") do |dep_erb|
      erb = ERB.new(dep_erb.read)
      erb.filename = dep_erb_path
      dep.print(erb.result)
    end
  end
end
