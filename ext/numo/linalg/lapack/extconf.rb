require 'mkmf'
require 'numo/narray'
require 'erb'

$LOAD_PATH.each do |x|
  if File.exist? File.join(x,'numo/numo/narray.h')
    $INCFLAGS = "-I#{x}/numo " + $INCFLAGS
    break
  end
end

srcs = %w(
lapack
lapack_s
lapack_d
lapack_c
lapack_z
)
$objs = srcs.collect{|i| i+".o"}

if !have_header('numo/narray.h')
  puts "
  Header numo/narray.h was not found. Give pathname as follows:
  % ruby extconf.rb --with-narray-include=narray_h_dir"
  exit(1)
end

if have_header("dlfcn.h")
  exit(1) unless have_library("dl")
  exit(1) unless have_func("dlopen")
elsif have_header("windows.h")
  exit(1) unless have_func("LoadLibrary")
end

dep_path = File.join(__dir__, "depend")
File.open(dep_path, "w") do |dep|
  dep_erb_path = File.join(__dir__, "depend.erb")
  File.open(dep_erb_path, "r") do |dep_erb|
    erb = ERB.new(dep_erb.read)
    erb.filename = dep_erb_path
    dep.print(erb.result)
  end
end

create_makefile('numo/linalg/lapack')
