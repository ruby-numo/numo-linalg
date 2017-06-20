require_relative "load_libs"

load_libs(use = ARGV[1])

a = Numo::DFloat.new(1500,1500).seq
b = Numo::DFloat.new(1500,1500).seq
n = 10

rep = ARGV[0].to_i
cap = ARGV[2]=="cap" ? Benchmark::CAPTION : ""
Benchmark.benchmark(cap,16) do |x|
  rep.times do
    x.report(use){ n.times{ Numo::Linalg::Blas.dgemm(a,b) } }
  end
end
