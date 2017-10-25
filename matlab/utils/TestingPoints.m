classdef TestingPoints < InterpolationPoints
    methods (Access = public)
        function Obj = TestingPoints(n_in_dims, n_pts, bounds, x_vals)
            Obj = Obj@InterpolationPoints(n_in_dims, n_pts, bounds);

            if (nargin > 3)
                Obj.recenterAndSave(x_vals);
            else
                test_pts = cell(Obj.n_in_dims, 1);
                for d_i = 1:Obj.n_in_dims
                    r_vec = rand(n_pts(d_i), 1);
                    cs_r_vec = cumsum(r_vec);
                    test_pts{d_i} = cs_r_vec;
                end

                % Each dimesion needs to be scaled and centered according to the original data dimesions, i.e., the
                % final data should be uniformly spread between the left and the right boundary: TODO
                Obj.recenterAndSave(test_pts);
            end
        end
    end
end
