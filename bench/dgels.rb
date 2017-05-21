require_relative "load_libs"

load_libs(use = ARGV[0])

a = Numo::DFloat.new(1000,1000).seq
b = Numo::DFloat.eye(1000)
n = 4

Benchmark.bm(16) do |x|
  x.report(use){ n.times{ Numo::Linalg::Lapack.dgels(a,b) } }
end
