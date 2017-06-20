require_relative "load_libs"

load_libs(use = ARGV[1])

a = Numo::DFloat.new(1000,1000).seq
b = Numo::DFloat.eye(1000)
n = 10

rep = ARGV[0].to_i
cap = ARGV[2]=="cap" ? Benchmark::CAPTION : ""
Benchmark.benchmark(cap,16) do |x|
  rep.times do
    x.report(use){ n.times{ Numo::Linalg::Lapack.dgels(a,b) } }
  end
end
