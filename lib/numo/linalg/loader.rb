require "fiddle"

module Numo
  module Linalg
    module Loader

      @@libs = nil

      module_function

      def libs
        @@libs
      end

      def load_mkl(path=nil,exc:true)
        a = case path
            when String; [path]
            when Array ;  path
            when NilClass; [MKL_LIBPATH,""].compact
            else raise ArgumentError,"invalid path"
            end
        a.each do |d|
          if d.empty?
            f_iomp5 = "libiomp5."+EXT
            f_mkl_core = "libmkl_core."+EXT
            f_intel_thread = "libmkl_intel_thread."+EXT
            f_intel_lp64 = "libmkl_intel_lp64."+EXT
          else
            e = d.sub(/\/mkl\//, "/")
            f_iomp5 = File.join(e,"libiomp5."+EXT)
            f_mkl_core = File.join(d,"libmkl_core."+EXT)
            f_intel_thread = File.join(d,"libmkl_intel_thread."+EXT)
            f_intel_lp64 = File.join(d,"libmkl_intel_lp64."+EXT)
          end
          begin
            Fiddle.dlopen(f_iomp5)
            Fiddle.dlopen(f_mkl_core)
            Fiddle.dlopen(f_intel_thread)
            Blas.dlopen(f_intel_lp64)
            Lapack.dlopen(f_intel_lp64)
            @@libs = [ f_iomp5, f_mkl_core, f_intel_thread, f_intel_lp64 ]
            if $DEBUG
              $stderr.puts "Numo::Linalg: use #{f_intel_lp64}"
            end
            return true
          rescue
          end
        end
        if exc
          raise RuntimeError, "cannot find MKL library"
        end
        false
      end

      def load_openblas(path=nil,exc:true)
        a = case path
            when String; [path]
            when Array ;  path
            when NilClass; [OPENBLAS_LIBPATH,""].compact
            else raise ArgumentError,"invalid path"
            end
        a.each do |d|
          %w[libopenblaso. libopenblasp. libopenblas.].each do |f|
            f = File.join(d,f) if !d.empty?
            f += EXT
            begin
              Blas.dlopen(f)
              Fiddle.dlopen(f)
              begin
                Lapack.dlopen("liblapacke."+EXT)
              rescue
                Lapack.dlopen(f)
              end
              @@libs = [ f ]
              if $DEBUG
                $stderr.puts "Numo::Linalg: use #{f}"
              end
              return true
            rescue
            end
          end
        end
        if exc
          raise RuntimeError, "cannot find OpenBLAS library"
        end
        false
      end

      def load_atlas(path=nil,exc:true)
        a = case path
            when String; [path]
            when Array ;  path
            when NilClass; [ATLAS_LIBPATH,""].compact
            else raise ArgumentError,"invalid path"
            end
        a.each do |d|
          %w[libtatlas. libatlas. libsatlas.].each do |f|
            f = File.join(d,f) if !d.empty?
            f += EXT
            begin
              Fiddle.dlopen(f)
              Blas.dlopen(f)
              l = LAPACK_LIBPATH
              e = "liblapacke."+EXT
              e = File.join(l,e) if l && !l.empty?
              Lapack.dlopen(e)
              @@libs = [ f, e ]
              if $DEBUG
                $stderr.puts "Numo::Linalg: use #{f}"
              end
              return true
            rescue
            end
          end
        end
        if exc
          raise RuntimeError, "cannot find ATLAS library"
        end
        false
      end

      def load_lapack(path=nil,exc:true)
        a = case path
            when String; [[path,path]]
            when Array ;  [path]
            when NilClass
              b = BLAS_LIBPATH || LAPACK_LIBPATH
              l = LAPACK_LIBPATH || BLAS_LIBPATH
              (b ? [[b,l]] : []) + ["",""]
            else raise ArgumentError,"invalid path"
            end
        a.each do |arg|
          b,l = *arg
          if b.empty?
            f_blas = "libblas."+EXT
            f_cblas = "libcblas."+EXT
            f_lapack = "liblapack."+EXT
            f_lapacke = "liblapacke."+EXT
          else
            f_blas = File.join(b,"libblas."+EXT)
            f_cblas = File.join(b,"libcblas."+EXT)
            f_lapack = File.join(l,"liblapack."+EXT)
            f_lapacke = File.join(l,"liblapacke."+EXT)
          end
          begin
            Fiddle.dlopen(f_blas)
            Fiddle.dlopen(f_lapack)
            Blas.dlopen(f_cblas)
            Lapack.dlopen(f_lapacke)
            @@libs = [ f_blas, f_lapack, f_cblas, f_lapacke ]
            if $DEBUG
              $stderr.puts "Numo::Linalg: use #{f_blas} and #{f_lapack}"
            end
            return true
          rescue
          end
        end
        if exc
          raise RuntimeError, "cannot find BLAS/LAPABK library"
        end
        false
      end

      def load_library
        case BACKEND
        when /mkl/i         ; return if load_mkl(exc:false)
        when /^openblas/i   ; return if load_openblas(exc:false)
        when /^atlas/i      ; return if load_atlas(exc:false)
        when /lapack|blas/i ; return if load_lapack(exc:false)
        else
          return if load_mkl(exc:false)
          return if load_openblas(exc:false)
          return if load_atlas(exc:false)
          return if load_lapack(exc:false)
        end
        raise RuntimeError, "cannot find backend library for Numo::Linalg"
      end

    end
  end
end
