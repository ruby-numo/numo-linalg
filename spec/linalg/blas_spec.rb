# frozen_string_literal: true

require 'spec_helper'

RSpec.describe Numo::Linalg::Blas do
  describe 'ger' do
    describe 'SFloat' do
      subject(:err) { (a - Numo::SFloat[[3, 4, 5], [6, 8, 10]]).sum.abs }

      let(:x) { Numo::SFloat[1, 2] }
      let(:y) { Numo::SFloat[3, 4, 5] }

      context 'when a is given' do
        let(:a) { Numo::SFloat.zeros(2, 3) }

        it 'calculates row and column vectors product' do
          Numo::Linalg::Blas.call(:ger, x, y, a)
          expect(err).to be < ERR_TOL
        end
      end

      context 'when a is not given' do
        let(:a) { Numo::Linalg::Blas.call(:ger, x, y) }

        it 'calculates row and column vectors product' do
          expect(err).to be < ERR_TOL
        end
      end
    end

    describe 'DFloat' do
      subject(:err) { (a - Numo::DFloat[[3, 4, 5], [6, 8, 10]]).sum.abs }

      let(:x) { Numo::DFloat[1, 2] }
      let(:y) { Numo::DFloat[3, 4, 5] }

      context 'when a is given' do
        let(:a) { Numo::DFloat.zeros(2, 3) }

        it 'calculates row and column vectors product' do
          Numo::Linalg::Blas.call(:ger, x, y, a)
          expect(err).to be < ERR_TOL
        end
      end

      context 'when a is not given' do
        let(:a) { Numo::Linalg::Blas.call(:ger, x, y) }

        it 'calculates row and column vectors product' do
          expect(err).to be < ERR_TOL
        end
      end
    end
  end

  describe 'syr2k' do
    describe 'SFloat' do
      subject(:err) { (c - Numo::SFloat[[50, 78, 106], [0, 90, 102], [0, 0, 98]]).sum.abs } # NOTE: syr2k returns only the upper triangular part of C.

      let(:a) { Numo::SFloat[[1, 2], [3, 4], [5, 6]] }
      let(:b) { Numo::SFloat[[9, 8], [7, 6], [5, 4]] }

      context 'when c is given' do
        let(:c) { Numo::SFloat.zeros(3, 3) }

        it 'performs a rank-2k matrix-matrix operation: C = A * B.T + B * A.T' do
          Numo::Linalg::Blas.call(:syr2k, a, b, c)
          expect(err).to be < ERR_TOL
        end
      end

      context 'when c is not given' do
        let(:c) { Numo::Linalg::Blas.call(:syr2k, a, b) }

        it 'performs a rank-2k matrix-matrix operation: A * B.T + B * A.T' do
          expect(err).to be < ERR_TOL
        end
      end
    end

    describe 'DFloat' do
      subject(:err) { (c - Numo::DFloat[[50, 78, 106], [0, 90, 102], [0, 0, 98]]).sum.abs }

      let(:a) { Numo::DFloat[[1, 2], [3, 4], [5, 6]] }
      let(:b) { Numo::DFloat[[9, 8], [7, 6], [5, 4]] }

      context 'when c is given' do
        let(:c) { Numo::DFloat.zeros(3, 3) }

        it 'performs a rank-2k matrix-matrix operation: C = A * B.T + B * A.T' do
          Numo::Linalg::Blas.call(:syr2k, a, b, c)
          expect(err).to be < ERR_TOL
        end
      end

      context 'when c is not given' do
        let(:c) { Numo::Linalg::Blas.call(:syr2k, a, b) }

        it 'performs a rank-2k matrix-matrix operation: A * B.T + B * A.T' do
          expect(err).to be < ERR_TOL
        end
      end
    end
  end
end
