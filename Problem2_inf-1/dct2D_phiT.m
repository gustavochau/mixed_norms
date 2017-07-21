function[v] = dct2D_phiT(b, Nrows, Ncols, factor)

    dummy = zeros( factor*Nrows*Ncols, 1 );
    v     = zeros( (factor*factor)*Nrows*Ncols, 1 );

    % Operate over columns

    n = Nrows;
    m = n*factor;
    f = (1:(m-1))';
    const = [n; .5 * (n +  sin(2.0*pi*f/factor) ./ (2.0*sin(pi*f/m)))];
    const = const .^ .5;

    n2 = 2*n;
    m = n * factor;
    y = zeros(4*m,1);

    for k=1:Ncols,
      
        y(2:2:n2) = b(1+(k-1)*Nrows:k*Nrows);
        z = fft(y);
        dummy(1+(k-1)*factor*Nrows:k*factor*Nrows) = real(z(1:m)) ./ const;
    end


    % Operate over rows

    n = Ncols;
    m = n*factor;
    f = (1:(m-1))';
    const = [n; .5 * (n +  sin(2.0*pi*f/factor) ./ (2.0*sin(pi*f/m)))];
    const = const .^ .5;

    n2 = 2*n;
    m = n * factor;
    y = zeros(4*m,1);

    for k=1:factor*Nrows,
      
        y(2:2:n2) = dummy(k:factor*Nrows:factor*Ncols*Nrows);
        z = fft(y);
        v(k:factor*Nrows:(factor*factor)*Ncols*Nrows) = real(z(1:m)) ./ const;
    end


%=============================================================
