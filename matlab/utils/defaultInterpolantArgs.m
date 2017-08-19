function args = defaultInterpolantArguments()
    args    = cell(4, 1);
    args{1} = 1; % Number of input dimensions. Input can only be a VECTOR
    args{2} = [-1; 1];
    args{3} = 4;
    args{4} = 'uniform';
end
