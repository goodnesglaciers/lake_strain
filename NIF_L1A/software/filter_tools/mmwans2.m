function C=mmwans2(A,B)

            Az = A ~= 0;
            Bz = B ~= 0;
            An = isnan(A);
            Bn = isnan(B);
            A(An) = 0;
            B(Bn) = 0;
            C = A * B;
            C(An * Bz | Az * Bn) = NaN;
