function[v] = dct2D_phi(b, Nrows, Ncols, factor)

    dummy = zeros( factor*Nrows*Ncols, 1 );
    v     = zeros( Nrows*Ncols, 1 );

    % Operate over columns

    m = factor*Nrows;
    n = Nrows;
    f = (1:(m-1))';
    const = [n; .5 * (n +  sin(2.0*pi*f/factor) ./ (2.0*sin(pi*f/m)))];
    const = const .^ .5;

    n2 = 2 * n;
    %b = b ./ const;
    z = zeros(4*m, 1);

    for k=1:factor*Ncols,
      z(1:m) = b(1+(k-1)*factor*Nrows:k*factor*Nrows) ./ const;
      y = fft(z);
      dummy( 1+(k-1)*Nrows:k*Nrows ) = real(y(2:2:n2));

    end


    % Operate over columns

    m = factor*Ncols;
    n = Ncols;
    f = (1:(m-1))';
    const = [n; .5 * (n +  sin(2.0*pi*f/factor) ./ (2.0*sin(pi*f/m)))];
    const = const .^ .5;

    n2 = 2 * n;
    %b = b ./ const;
    z = zeros(4*m, 1);


    for k=1:Nrows,

      z(1:m) = dummy(k:Nrows:factor*Nrows*Ncols) ./ const;
      y = fft(z);
      v( k:Nrows:Nrows*Ncols ) = real(y(2:2:n2));

    end

%=============================================================
%=============================================================

