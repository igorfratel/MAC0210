# Allan Rocha e Igor Fratel
# n. 9761614 | n. 9793565
# Interpolação polinomial por partes bivariada
1;

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
# fmatrix = matriz com os valores da função f em cada ponto da malha
# dxfmatrix = matriz com os valores da derivada parcial em x da f em cada ponto da malha
# dyfmatrix = matriz com os valores da derivada parcial em y da f em cada ponto da malha
# dxyfmatrix = matriz com os valores da derivada segunda parcial em x e y (ou vice versa) em cada ponto da malha
# typeintpol = tipo de interpolação: bilinear ou bicúbo
#
# A funçao retorna:
# Uma matriz com os coeficientes para cada ponto da malha para que se possa fazer a função interpoladora
function coefficients = constroiv(nx, ny, ax, bx, ay, by, fmatrix, dxfmatrix, dyfmatrix, dxyfmatrix, typeintpol)
    # Configurações iniciais:
    hx = (bx - ax)/nx;
    hy = (by - ay)/ny;
    fmatrix = transpose(fmatrix);
    dxfmatrix = transpose(dxfmatrix);
    dyfmatrix = transpose(dyfmatrix);
    dxyfmatrix = transpose(dxyfmatrix);

    if (strcmp(typeintpol, "linear"))
        coefficients = zeros(2*nx, 2*ny);
        i = 1;
        ic = 1;
        while (i <= nx)
            j = 1;
            jc = 1;
            while (j <= ny)
                # fmatrix(i, j) ==> f(x_i, y_j)
                # fmatrix(i + 1, j) ==> f(x_i+1, y_j)
                # ...

                # Coeficiente elementar
                a00 = fmatrix(i, j);

                # Preenchendo a matriz de coeficientes
                coefficients(ic, jc) = a00; # a00
                a01 = fmatrix(i, j+1) - a00; # Pseudo-elementar
                coefficients(ic, jc+1) = a01; # a01
                a10 = fmatrix(i+1, j) - a00; # Pseudo-elementar
                coefficients(ic+1, jc) = a10; # a10
                coefficients(ic+1, jc+1) = fmatrix(i+1, j+1) - a00 - a10 - a01; # a11
                j += 1;
                jc += 2;
            endwhile
            i += 1;
            ic += 2;
        endwhile
    endif

    if (strcmp(typeintpol, "cubico"))
        coefficients = zeros(4*nx, 4*ny);
        i = 1;
        ic = 1;
        while (i <= nx)
            j = 1;
            jc = 1;
            while (j <= nx)
                # dxfmatrix(i, j) ==> dxf(x_i, y_j) - Derivada parcial no x
                # dyfmatrix(i + 1, j) ==> dyf(x_i+1, y_j) - Derivada parcial no y
                # ...

                # Coeficientes elementares
                a00 = fmatrix(i, j);
                a01 = hy*dyfmatrix(i, j);
                a10 = hx*dxfmatrix(i, j);
                a11 = hx*hy*dxyfmatrix(i, j);

                # Preenchendo os coeficientes
                coefficients(ic, jc) = a00; # a00
                coefficients(ic, jc+1) = a01; # a01
                a02 = 3*fmatrix(i, j+1) - hy*dyfmatrix(i, j+1) -2*hy*dyfmatrix(i, j) - 3*a00; # Pseudo-elementar
                coefficients(ic, jc+2) = a02; # a02
                a03 = fmatrix(i, j+1) - a00 - a01 - a02; # Pseudo-elementar
                coefficients(ic, jc+3) = a03; # a03

                coefficients(ic+1, jc) = a10; # a10
                coeficientes(ic+1, jc+1) = a11; # a11
                a12 = 3*hx*dxfmatrix(i, j+1) - hx*hy*dxyfmatrix(i, j+1) - 3*a10 - 2*a11; # Pseudo-elementar
                coefficients(ic+1, jc+2) = a12; # a12
                a13 = hx*dxfmatrix(i, j+1) - a10 - a11 - a12; # Pseudo-elementar
                coefficients(ic+1, jc+3) = a13; # a13

                a20 = 3*fmatrix(i+1, j) - hx*dxfmatrix(i+1, j) - 2*hx*dxfmatrix(i, j) - 3*a00; # Pseudo-elementar
                coefficients(ic+2, jc) = a20; # a20
                a21 = 3*hy*dyfmatrix(i+1, j) - hx*hy*dxyfmatrix(i+1, j) - 3*a01 - 2*a11; # Pseudo-elementar
                coefficients(ic+2, jc+1) = a21; # a21

                a30 = fmatrix(i+1, j) - a00 - a10 - a20; # Pseudo-elementar
                coefficients(ic+3, jc) = a30; # a30
                a31 = hy*dyfmatrix(i+1, j) - a01 - a11 - a21; # Pseudo-elementar
                coefficients(ic+3, jc+1) = a31; # a31
                
                # Somas auxiliares constantes
                c = a00 + a01 + a02 + a03 + a10 + a11 + a12 + a13 + a20 + a21 + a30 + a31;
                d = a10 + a11 + a12 + a13 + 2*(a20 + a21) + 3*(a30 + a31);
                e = a01 + a11 + a21 + a31 + 2*(a02 + a12) + 3*(a03 + a13);
                f = a11 + 2*a21 + 3*a31 + 2*a12 + 3*a13;

                matrice_inverse = [9, -3, -3, 1; -6, 2, 3, -1; -6, 3, 2, -1; 4, -2, -2, 1];
                matrice_functions = [fmatrix(i+1, j+1) - c; hx*dxfmatrix(i+1, j+1) - d; hy*dyfmatrix(i+1, j+1) - e; hx*hy*dxyfmatrix(i+1, j+1) - f];

                inter_coef = matrice_inverse*matrice_functions;

                coefficients(ic+2, jc+2) = inter_coef(1); # a22
                coefficients(ic+2, jc+3) = inter_coef(2); # a23
                coefficients(ic+3, jc+2) = inter_coef(3); # a32
                coefficients(ic+3, jc+3) = inter_coef(4); # a33
                j += 1;
                jc += 4;
            endwhile
            i += 1;
            ic += 4;
        endwhile
    endif

endfunction

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
# coefficients = matriz com os coeficientes necessários para fazer a função interpoladora em cada ponto da malha
# typeintpol = tipo de interpolação: bilinear ou bicúbo
#
# A funçao retorna:
# O valor da função interpoladora no ponto (x, y)
function s = avaliav(nx, ny, ax, bx, ay, by, x, y, coefficients, typeintpol)
    # Configurações iniciais
    hx = (bx - ax)/nx;
    hy = (by - ay)/ny;

    # Índices
    i = floor((x - ax)/hx);
    j = floor((y - ay)/hy);

    # Valores de x_i e y_j
    xi = ax + i*hx;
    yj = ay + j*hy;

    # Correção nos índices
    if (i != nx) i++;
    else xi -= hx;
    endif
    if (j != ny) j++;
    else yj -= hy;
    endif

    # Alfa e beta
    alfa = (x - xi)/hx;
    beta = (y - yj)/hy;

    if (strcmp(typeintpol, "linear"))
        # Índices deslocados
        i = 2*i - 1;
        j = 2*j - 1;

        s = coefficients(i, j) + coefficients(i+1, j)*alfa + coefficients(i, j+1)*beta + coefficients(i+1, j+1)*alfa*beta;
    endif

    if (strcmp(typeintpol, "cubico"))
        # Índices deslocados
        i = 4*i - 3;
        j = 4*j - 3;

        vector_x = [1, alfa, alfa^2, alfa^3];
        vector_y = [1; beta; beta^2; beta^3];
        coeff = [coefficients(i, j), coefficients(i, j+1), coefficients(i, j+2), coefficients(i, j+3); coefficients(i+1, j), coefficients(i+1, j+1), coefficients(i+1, j+2), coefficients(i+1, j+3); coefficients(i+2, j), coefficients(i+2, j+1), coefficients(i+2, j+2), coefficients(i+2, j+3); coefficients(i+3, j), coefficients(i+3, j+1), coefficients(i+3, j+2), coefficients(i+3, j+3)];

        s = vector_x*coeff*vector_y;
    endif

endfunction