# frozen_string_literal: true

require 'spec_helper'

RSpec.describe Numo::Linalg do
  describe 'lu_fact' do
    let(:m) { 5 }
    let(:n) { 3 }
    let(:mat_a) { rand_rect_real_mat(m, n) }
    let(:mat_b) { rand_rect_complex_mat(m, n) }

    def permutation_mat(ipiv)
      mat_p = Numo::DFloat.eye(m)
      (ipiv - 1).to_a.each_with_index.select do |a, b|
        unless a - b == 0
          tmp = mat_p[b, true].dup
          mat_p[b, true] = mat_p[a, true]
          mat_p[a, true] = tmp
        end
      end
      mat_p
    end

    it 'raises ArgumentError given a invalid mode option' do
      expect { described_class.lu_fact(Numo::DFloat.new(2, 4).rand, mode: 'foo') }.to raise_error(ArgumentError)
    end

    it 'calculates the LU factorization of a rectangular real matrix' do
      lu, ipiv = described_class.lu_fact(mat_a)
      mat_l = lu.tril
      mat_l[mat_l.diag_indices(0)] = 1.0
      mat_u = lu.triu[0...n, 0...n]
      mat_p = permutation_mat(ipiv)
      expect((mat_p.dot(mat_a) - mat_l.dot(mat_u)).abs.max).to be < ERR_TOL
    end

    it 'calculates the LU factorization of a rectangular complex matrix' do
      lu, ipiv = described_class.lu_fact(mat_b)
      mat_l = lu.tril
      mat_l[mat_l.diag_indices(0)] = 1.0
      mat_u = lu.triu[0...n, 0...n]
      mat_p = permutation_mat(ipiv)
      expect((mat_p.dot(mat_b) - mat_l.dot(mat_u)).abs.max).to be < ERR_TOL
    end
  end
end
