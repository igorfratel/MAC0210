# Allan Rocha e Igor Fratel
# n. 9761614 | n. 9793565
# Interpolação e diferenciação numérica
1;

# VARIÁVEIS GLOBAIS PARA O AVALIAV

global NX;
global NY;
global AX;
global BX;
global AY;
global BY;
global COEFFICIENTS;

# FUNÇÕES USADAS NOS TESTES

# sen(xy)
function z = sen(x, y)
    z = sin(x*y);
endfunction

# y.sen(x)
function z = ysen(x, y)
    z = y*sin(x);
endfunction

# x + y
function z = f1(x, y)
    z = x + y;
endfunction

# x² + y²
function z = f2(x, y)
    z = x^2 + y^2;
endfunction

# Derivada parcial em x de y.sen(x)
function w = ysenx(x, y)
    w = cos(x)*y;
endfunction

# Derivada parcial em y de y.sen(x)
function u = yseny(x, y)
    u = sin(x);
endfunction

# Derivada parcial em x e y de y.sen(x)
function v = ysenxy(x, y)
    v = cos(x);
endfunction

# FUNÇÕES AUXILIARES PARA IMPRESSÃO/EXECUÇÃO DOS TESTES, ETC

function numerical = imagem2(nx, ny, ax, bx, ay, by, foto, n)
    p = 2;
    while (p <= columns(foto))
        foto(:,[p]) = [];
        p += 1;
    endwhile
    foto = transpose(foto);
    p = 2;
    while (p <= columns(foto))
        foto(:,[p]) = [];
        p += 1;
    endwhile
    foto = transpose(foto);

    # foto = double(foto);

    coef = constroiv(nx/(2*n), ny/(2*n), ax, bx, ay, by, foto, "cubico");

    numerical = zeros(nx + 1, ny + 1);

    hx = (bx - ax)/nx;
    hy = (by - ay)/ny;

    i = 0;
    while (i <= nx)
        j = 0;
        while (j <= ny)
            xi = ax + i*hx;
            yj = ay + j*hy;
            numerical(i+1, j+1) = avaliav(xi, yj);
            j += 1;
        endwhile
        i += 1;
    endwhile

    numerical = transpose(numerical);

endfunction

function imagem(nx, ny, ax, bx, ay, by, f, n)
    m = fm(nx, ny, ax, bx, ay, by, f);
    m1 = fm(nx/(2*n), ny/(2*n), ax, bx, ay, by, f);
    coef = constroiv(nx/(2*n), ny/(2*n), ax, bx, ay, by, m1, "cubico");

    numerical = zeros(nx + 1, ny + 1);

    hx = (bx - ax)/nx;
    hy = (by - ay)/ny;

    i = 0;
    while (i <= nx)
        j = 0;
        while (j <= ny)
            xi = ax + i*hx;
            yj = ay + j*hy;
            numerical(i+1, j+1) = avaliav(xi, yj);
            j += 1;
        endwhile
        i += 1;
    endwhile

    numerical = transpose(numerical);


    disp(m);
    disp("");
    disp(m1);
    disp("");
    disp(numerical);
    disp("");

endfunction

function me = interpola(nx, ny, ax, bx, ay, by, f)
    m = fm(nx, ny, ax, bx, ay, by, f);
    coef = constroiv(nx, ny, ax, bx, ay, by, m, "cubico");
    analytical = fm(nx*2, ny*2, ax, bx, ay, by, f);
    numerical = zeros(nx*2 + 1, ny*2 + 1);

    h = (bx - ax)/(nx*2);
    k = (by - ay)/(ny*2);

    i = 0;
    while (i <= nx*2)
        j = 0;
        while (j <= ny*2)
            xi = ax + i*h;
            yj = ay + j*k;
            numerical(i+1, j+1) = avaliav_(nx, ny, ax, bx, ay, by, xi, yj, coef, "cubico");
            j += 1;
        endwhile
        i += 1;
    endwhile

    numerical = transpose(numerical);

    me = abs(analytical - numerical);

    disp("Pontos a interpolar:")
    disp("");
    disp(m);
    disp("");
    disp("Pontos analíticos (f)");
    disp("");
    disp(analytical);
    disp("");
    disp("Pontos numéricos (v)");
    disp("");
    disp(numerical);
    disp("");
    disp("| f - v |");
    disp("");
    disp(me);
    disp("");

endfunction

# A função recebe:
# error1 = matriz de erros
# erro2 = matriz de erros
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
#
# A função retorna:
# O fator de erro entre as duas matrizes de erros
function e = lenerror(error1, error2, nx, ny)
    n = (nx+1)*(ny+1);
    mediaI = 0;
    mediaII = 0;

    i = 1;
    while (i <= nx + 1)
        j = 1;
        while (j <= ny + 1)
                mediaI += error1(i, j);
                mediaII += error2(i, j);
            j += 1;
        endwhile
        i += 1;
    endwhile

    mediaI /= n;
    mediaII /= n;

    e = mediaI/mediaII;

    disp(e);

endfunction

# A função recebe:
# typed = tipo de derivada parcial
# fmatrix = matriz com os valores da função f em cada ponto da malha
# mh = x/y máximo
#
# A função retorna:
# Uma matriz com o erro da derivada parcial real e aproximada de sen(xy)
function ematrix = error(typed, fmatrix, mh)
    fxr = [];
    fxa = [];
    if (strcmp(typed, "x"))
        fxr = fm(5, 5, 0, mh, 0, mh, @ysenx);
        fxa = aproxdf(5, 5, 0, mh, 0, mh, fmatrix, "x");
    else
        if (strcmp(typed, "y"))
            fxr = fm(5, 5, 0, mh, 0, mh, @yseny);
            fxa = aproxdf(5, 5, 0, mh, 0, mh, fmatrix, "y");
        else
            if (strcmp(typed, "xy"))
                fxr = fm(5, 5, 0, mh, 0, mh, @ysenxy);
                fy = aproxdf(5, 5, 0, mh, 0, mh, fmatrix, "y");
                fxa = aproxdf(5, 5, 0, mh, 0, mh, fy, "xy");
            endif
        endif
    endif

    if (strcmp(typed, "xy"))
        disp("Malha de derivadas parciais de segunda ordem mistas reais em x e y");
        disp("");
        disp(fxr);
        disp("");
        disp("Malha de derivadas parciais de segunda ordem mistas aproximadas em x e y");
        disp("");
        disp(fxa);
        disp("");
    else
        disp(cstrcat("Malha de derivadas parciais reais em ", typed));
        disp("");
        disp(fxr);
        disp("");
        disp(cstrcat("Malha de derivadas parciais aproximadas em ", typed));
        disp("");
        disp(fxa);
        disp("");
    endif

    ematrix = abs(fxr - fxa);
endfunction

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
#
# A função retorna:
# Uma matriz representando a função 'f' nos pontos da malha em questão
function matrix = fm(nx, ny, ax, bx, ay, by, f)
    # Configurações iniciais:
    h = (bx - ax)/nx;
    k = (by - ay)/ny;
    matrix = zeros(nx + 1, ny + 1);

    i = 0;
    while (i <= nx)
        j = 0;
        while (j <= ny)
            xi = ax + i*h;
            yj = ay + j*k;
            matrix(i+1, j+1) = f(xi, yj);
            j += 1;
        endwhile
        i += 1;
    endwhile

    matrix = transpose(matrix);

endfunction

### ### ### ### ### ### ### ###

### FUNÇÕES RELATIVAS AO EP ###

### ### ### ### ### ### ### ###

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
# fmatrix = matriz com os valores da função f em cada ponto da malha ou da derivada parcial da f em relação a x em cada ponto da malha
# typed = tipo de derivada: "x", "y" ou "xy"
#
# A função retorna:
# Uma matriz com os valores da função fx, fy ou fxy em cada ponto da malha.
function dfmatrix = aproxdf(nx, ny, ax, bx, ay, by, fmatrix, typed)
    # Configurações iniciais:
    h = (bx - ax)/nx;
    k = (by - ay)/ny;
    fmatrix = transpose(fmatrix);
    dfmatrix = zeros(nx + 1, ny + 1);

    if (strcmp(typed, "x"))
        i = 1;
        while (i <= nx + 1)
            j = 1;
            while (j <= ny + 1)
                # fmatrix(i, j) ==> f(x_i, y_j)
                # fmatrix(i + 1, j) ==> f(x_i+1, y_j)
                # ...

                if (i == 1)
                    dfmatrix(i, j) = (4*fmatrix(i+1, j) - fmatrix(i+2, j) - 3*fmatrix(i, j))/(2*h);
                else
                    if (i == nx + 1)
                        dfmatrix(i, j) = (-4*fmatrix(i-1, j) + fmatrix(i-2, j) + 3*fmatrix(i, j))/(2*h);
                    else
                        dfmatrix(i, j) = (fmatrix(i+1, j) - fmatrix(i-1, j))/(2*h);
                    endif
                endif
                j += 1;
            endwhile
            i += 1;
        endwhile
    endif

    if (strcmp(typed, "y"))
        i = 1;
        while (i <= nx + 1)
            j = 1;
            while (j <= ny + 1)
                # fmatrix(i, j) ==> f(x_i, y_j)
                # fmatrix(i + 1, j) ==> f(x_i+1, y_j)
                # ...

                if (j == 1)
                    dfmatrix(i, j) = (4*fmatrix(i, j+1) - fmatrix(i, j+2) - 3*fmatrix(i, j))/(2*k);
                else
                    if (j == ny + 1)
                        dfmatrix(i, j) = (-4*fmatrix(i, j-1) + fmatrix(i, j-2) + 3*fmatrix(i, j))/(2*k);
                    else
                        dfmatrix(i, j) = (fmatrix(i, j+1) - fmatrix(i, j-1))/(2*k);
                    endif
                endif
                j += 1;
            endwhile
            i += 1;
        endwhile
    endif

    if (strcmp(typed, "xy"))
        i = 1;
        while (i <= nx + 1)
            j = 1;
            while (j <= ny + 1)
                # fmatrix(i, j) ==> f(x_i, y_j)
                # fmatrix(i + 1, j) ==> f(x_i+1, y_j)
                # ...

                if (i == 1)
                    dfmatrix(i, j) = (4*fmatrix(i+1, j) - fmatrix(i+2, j) - 3*fmatrix(i, j))/(2*h);
                else
                    if (i == nx + 1)
                        dfmatrix(i, j) = (-4*fmatrix(i-1, j) + fmatrix(i-2, j) + 3*fmatrix(i, j))/(2*h);
                    else
                        dfmatrix(i, j) = (fmatrix(i+1, j) - fmatrix(i-1, j))/(2*h);
                    endif
                endif
                j += 1;
            endwhile
            i += 1;
        endwhile
    endif

    dfmatrix = transpose(dfmatrix);

endfunction

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
# fmatrix = matriz com os valores da função f em cada ponto da malha
# typeintpol = tipo de interpolação: (bi)"linear" ou (bi)"cubico"
#
# A função retorna:
# Uma matriz com os coeficientes para cada ponto da malha para que se possa fazer a função interpoladora
function coefficients = constroiv(nx, ny, ax, bx, ay, by, fmatrix, typeintpol)
    # Setando variáveis globais
    global NX;
    global NY;
    global AX;
    global BX;
    global AY;
    global BY;
    global COEFFICIENTS;
    NX = nx;
    NY = ny;
    AX = ax;
    BX = bx;
    AY = ay;
    BY = by;
    # Configurações iniciais:
    hx = (bx - ax)/nx;
    hy = (by - ay)/ny;
    fmatrix = transpose(fmatrix);
    dxfmatrix = transpose(aproxdf(nx, ny, ax, bx, ay, by, fmatrix, "x"));
    dyfmatrix = transpose(aproxdf(nx, ny, ax, bx, ay, by, fmatrix, "y"));
    dxyfmatrix = transpose(aproxdf(nx, ny, ax, bx, ay, by, transpose(dyfmatrix), "xy"));

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
                a03 = 2*(fmatrix(i, j) - fmatrix(i, j+1)) + hy*(dyfmatrix(i, j) + dyfmatrix(i, j+1)); # Pseudo-elementar
                coefficients(ic, jc+3) = a03; # a03

                coefficients(ic+1, jc) = a10; # a10
                coeficientes(ic+1, jc+1) = a11; # a11
                a12 = 3*hx*dxfmatrix(i, j+1) - hx*hy*dxyfmatrix(i, j+1) - 3*a10 - 2*a11; # Pseudo-elementar
                coefficients(ic+1, jc+2) = a12; # a12
                a13 = 2*hx*(dxfmatrix(i, j) - dxfmatrix(i, j+1)) + hx*hy*(dxyfmatrix(i, j) + dxyfmatrix(i, j+1)); # Pseudo-elementar
                coefficients(ic+1, jc+3) = a13; # a13

                a20 = 3*fmatrix(i+1, j) - hx*dxfmatrix(i+1, j) - 2*hx*dxfmatrix(i, j) - 3*a00; # Pseudo-elementar
                coefficients(ic+2, jc) = a20; # a20
                a21 = 3*hy*dyfmatrix(i+1, j) - hx*hy*dxyfmatrix(i+1, j) - 3*a01 - 2*a11; # Pseudo-elementar
                coefficients(ic+2, jc+1) = a21; # a21

                a30 = 2*(fmatrix(i, j) - fmatrix(i+1, j)) + hx*(dxfmatrix(i, j) + dxfmatrix(i+1, j)); # Pseudo-elementar
                coefficients(ic+3, jc) = a30; # a30
                a31 = 2*hy*(dyfmatrix(i, j) - dyfmatrix(i+1, j)) + hx*hy*(dxyfmatrix(i, j) + dxyfmatrix(i+1, j)); # Pseudo-elementar
                coefficients(ic+3, jc+1) = a31; # a31
                
                # Somas auxiliares constantes
                c = a00 + a01 + a02 + a03 + a10 + a11 + a12 + a13 + a20 + a21 + a30 + a31;
                d = hx*(a10 + a11 + a12 + a13 + 2*a20 + 2*a21 + 3*a30 + 3*a31);
                e = hy*(a01 + a11 + a21 + a31 + 2*a02 + 2*a12 + 3*a03 + 3*a13);
                f = hx*hy*(a11 + 2*a21 + 2*a12 + 3*a31 + 9*a13);

                matrice_inverse = [9, -3, -3, 1; -6, 2, 3, -1; -6, 3, 2, -1; 4, -2, -2, 1];
                matrice_functions = [fmatrix(i+1, j+1) - c; (dxfmatrix(i+1, j+1) - d)/hx; (dyfmatrix(i+1, j+1) - e)/hy; (dxyfmatrix(i+1, j+1) - f)/(hx*hy)];

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

    COEFFICIENTS = coefficients;

endfunction

# A função recebe:
# nx = número de intervalos no eixo x
# ny = número de intervalos no eixo y
# ax = x mínimo
# bx = x máximo
# ay = y mínimo
# by = y máximo
# x = x a ser avaliado
# y = y a ser avaliado
# coefficients = matriz com os coeficientes necessários para fazer a função interpoladora em cada ponto da malha
# typeintpol = tipo de interpolação: bilinear ou bicúbica
#
# A funçao retorna:
# O valor da função interpoladora no ponto (x, y)
function s = avaliav_(nx, ny, ax, bx, ay, by, x, y, coefficients, typeintpol)
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

# A função recebe:
# x = x a ser avaliado
# y = y a ser avaliado
#
# A funçao retorna:
# O valor da função interpoladora no ponto (x, y)
function s = avaliav(x, y)
    global NX;
    global NY;
    global AX;
    global BX;
    global AY;
    global BY;
    global COEFFICIENTS;
    s = avaliav_(NX, NY, AX, BX, AY, BY, x, y, COEFFICIENTS, "cubico");
endfunction

function main(h)
    disp("***TESTE I - APROXIMAÇÃO DAS DERIVADAS***");
    disp("___Função utilizada: y.sen(x)___")
    disp("");

    h *= 5;

    fmatrix = fm(5, 5, 0, h, 0, h, @ysen);
    disp(cstrcat("(malha com h = ", num2str(h/5),")"));
    disp("");

    disp("DERIVADA PARCIAL EM X");
    disp("");
    errorsenxI = error("x", fmatrix, h);
    disp("Malha de erro entre derivadas parciais (em x) reais e aproximadas");
    disp("");
    disp(errorsenxI);
    disp("");

    disp("DERIVADA PARCIAL EM Y")
    disp("");
    errorsenyI = error("y", fmatrix, h);
    disp("Malha de erro entre derivadas parciais (em y) reais e aproximadas");
    disp("");
    disp(errorsenyI);
    disp("");

    disp("DERIVADA PARCIAL DE SEGUNDA ORDEM MISTAS EM X E Y")
    disp("");
    errorsenxyI = error("xy", fmatrix, h);
    disp("Malha de erro entre derivadas parciais de segunda ordem mistas (em x e y) reais e aproximadas");
    disp("");
    disp(errorsenxyI);
    disp("");

    fmatrix = fm(5, 5, 0, h/2, 0, h/2, @ysen);
    disp(cstrcat("(malha com h = ", num2str(h/5*(1/2)),")"));
    disp("");

    disp("DERIVADA PARCIAL EM X")
    disp("");
    errorsenxII = error("x", fmatrix, h/2);
    disp("Malha de erro entre derivadas parciais (em x) reais e aproximadas");
    disp("");
    disp(errorsenxII);
    disp("");

    disp("DERIVADA PARCIAL EM Y")
    disp("");
    errorsenyII = error("y", fmatrix, h/2);
    disp("Malha de erro entre derivadas parciais (em y) reais e aproximadas");
    disp("");
    disp(errorsenyII);
    disp("");

    disp("DERIVADA PARCIAL DE SEGUNDA ORDEM MISTAS EM X E Y")
    disp("");
    errorsenxyII = error("xy", fmatrix, h/2);
    disp("Malha de erro entre derivadas parciais de segunda ordem mistas (em x e y) reais e aproximadas");
    disp("");
    disp(errorsenxyII);
    disp("");

    disp("FATORES DE ERROS");
    disp("(O quanto o erro diminui quando h é dividido por dois)");
    disp("");
    disp("Fator de erro da derivada parcial em x:");
    fe1 = lenerror(errorsenxI, errorsenxII, 5, 5);
    disp("");
    disp("Fator de erro da derivada parcial em y:");
    fe2 = lenerror(errorsenyI, errorsenyII, 5, 5);
    disp("");
    disp("Fator de erro da derivada parcial de segunda ordem mistas em x e y:")
    fe3 = lenerror(errorsenxyI, errorsenxyII, 5, 5);
    disp("")
    disp("Fator de erro das derivadas parciais:")
    disp((fe1+fe2+fe3)/3);
    disp("")

    disp("***TESTE II - INTERPOLAÇÃO***");
    disp("");
    disp("--> FUNÇÃO - X + Y")
    disp("");
    interpola(2, 2, 0, h, 0, h, @f1);
    disp("--> FUNÇÃO - X² + X²")
    disp("");
    interpola(2, 2, 0, h, 0, h, @f2);
    disp("--> FUNÇÃO - SEN(XY)")
    disp("");
    interpola(2, 2, 0, h, 0, h, @sen);

    disp("***TESTE III - COMPRESSÃO DE IMAGEM***");
    disp("");

    # CONVERTENDO IMAGEM EM MATRIZ
    matrix = imread("mario.png");

    R = double(matrix(:,:,1));
    G = double(matrix(:,:,2));
    B = double(matrix(:,:,3));

    # ELIMINANDO COLUNA E LINHA DO R
    R(:,[16]) = [];
    R = transpose(R);
    R(:,[16]) = [];
    R = transpose(R);

    # ELIMINANDO COLUNA E LINHA DO G
    G(:,[16]) = [];
    G = transpose(G);
    G(:,[16]) = [];
    G = transpose(G);

    # ELIMINANDO COLUNA E LINHA DO R
    B(:,[16]) = [];
    B = transpose(B);
    B(:,[16]) = [];
    B = transpose(B);

    # COMPRIMINDO R
    R2 = imagem2(14, 14, 0, 15, 0, 15, R, 1);

    # COMPRIMINDO G
    G2 = imagem2(14, 14, 0, 15, 0, 15, G, 1);

    # COMPRIMINDO B
    B2 = imagem2(14, 14, 0, 15, 0, 15, B, 1);

    # DESCOMPRIMINDO
    foto2 = cat(3, R2, G2, B2);
    foto2 = uint8(foto2);
    disp("Mario");
    imshow(foto2);

endfunction

main(0.6)