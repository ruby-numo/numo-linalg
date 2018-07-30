# frozen_string_literal: true

require 'spec_helper'

RSpec.describe Numo::Linalg do
  describe 'solve' do
    let(:m) { 5 }
    let(:n) { 3 }
    let(:mat_a) { rand_square_real_mat(m) }
    let(:mat_b) { rand_rect_real_mat(m, n) }
    let(:vec_b) { rand_real_vec(m) }
    let(:mat_c) { rand_square_complex_mat(m) }
    let(:mat_d) { rand_rect_complex_mat(m, n) }
    let(:vec_d) { rand_complex_vec(m) }

    it 'raises ArgumentError given a invalid driver option' do
      expect { described_class.solve(mat_a, vec_b, driver: 'foo') }.to raise_error(ArgumentError)
    end

    it 'raises ShapeError when the number of rows of matrices are different' do
      expect { described_class.solve(mat_a, mat_b.transpose) }.to raise_error(Numo::NArray::ShapeError)
    end

    it 'raises ShapeError given a rectangular matrix as matrix A' do
      expect { described_class.solve(mat_b, mat_a) }.to raise_error(Numo::NArray::ShapeError)
    end

    it 'solves the linear equation A x = b with a square real matrix A' do
      vec_x = described_class.solve(mat_a, vec_b)
      expect((mat_a.dot(vec_x) - vec_b).abs.max).to be < ERR_TOL
      mat_x = described_class.solve(mat_a, mat_b)
      expect((mat_a.dot(mat_x) - mat_b).abs.max).to be < ERR_TOL
      vec_x = described_class.solve(mat_a, vec_d)
      expect((mat_a.dot(vec_x) - vec_d).abs.max).to be < ERR_TOL
      mat_x = described_class.solve(mat_a, mat_d)
      expect((mat_a.dot(mat_x) - mat_d).abs.max).to be < ERR_TOL
    end

    it 'solves the linear equation A x = b with a square complex matrix A' do
      vec_x = described_class.solve(mat_c, vec_b)
      expect((mat_c.dot(vec_x) - vec_b).abs.max).to be < ERR_TOL
      mat_x = described_class.solve(mat_c, mat_b)
      expect((mat_c.dot(mat_x) - mat_b).abs.max).to be < ERR_TOL
      vec_x = described_class.solve(mat_c, vec_d)
      expect((mat_c.dot(vec_x) - vec_d).abs.max).to be < ERR_TOL
      mat_x = described_class.solve(mat_c, mat_d)
      expect((mat_c.dot(mat_x) - mat_d).abs.max).to be < ERR_TOL
    end
  end
end
