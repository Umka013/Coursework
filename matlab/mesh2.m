clc;
clear;
close all;

elemSize = [1 1];
a = elemSize(1);
b = elemSize(2);

MeshNodes = [
    0 0;
    1 0;
    2 0;
    3 0;
    4 0;
    5 0;
    6 0;
    7 0;
    8 0;

    0 1;
    1 1;
    2 1;
    3 1;
    4 1;
    5 1;
    6 1;
    7 1;
    8 1;

    0 2;
    1 2;
    2 2;
    3 2;
    4 2;
    5 2;
    6 2;
    7 2;
    8 2;

    0 3;
    1 3;
    2 3;
    3 3;
    4 3;
    5 3;
    6 3;
    7 3;
    8 3;

    0 4;
    1 4;
    2 4;
    3 4;
    4 4;
    5 4;
    6 4;
    7 4;
    8 4;

    0 5;
    1 5;
    2 5;
    3 5;
    4 5;
    5 5;
    6 5;
    7 5;
    8 5;

    0 6;
    1 6;
    2 6;
    3 6;
    4 6;
    5 6;
    6 6;
    7 6;
    8 6;

    0 7;
    1 7;
    2 7;
    3 7;
    4 7;
    5 7;
    6 7;
    7 7;
    8 7;

    0 8;
    1 8;
    2 8;
    3 8;
    4 8;
    5 8;
    6 8;
    7 8;
    8 8;
];

MeshElems = [
    1   2   11  10;
    2   3   12  11;
    3   4   13  12;
    4   5   14  13;
    5   6   15  14;
    6   7   16  15;
    7   8   17  16;
    8   9   18  17;

    10  11  20  19;
    11  12  21  20;
    12  13  22  21;
    13  14  23  22;
    14  15  24  23;
    15  16  25  24;
    16  17  26  25;
    17  18  27  26;

    19  20  29  28;
    20  21  30  29;
    21  22  31  30;
    22  23  32  31;
    23  24  33  32;
    24  25  34  33;
    25  26  35  34;
    26  27  36  35;

    28  29  38  37;
    29  30  39  38;
    30  31  40  39;
    31  32  41  40;
    32  33  42  41;
    33  34  43  42;
    34  35  44  43;
    35  36  45  44;

    37  38  47  46;
    38  39  48  47;
    39  40  49  48;
    40  41  50  49;
    41  42  51  50;
    42  43  52  51;
    43  44  53  52;
    44  45  54  53;

    46  47  56  55;
    47  48  57  56;
    48  49  58  57;
    49  50  59  58;
    50  51  60  59;
    51  52  61  60;
    52  53  62  61;
    53  54  63  62;

    55  56  65  64;
    56  57  66  65;
    57  58  67  66;
    58  59  68  67;
    59  60  69  68;
    60  61  70  69;
    61  62  71  70;
    62  63  72  71;

    64  65  74  73;
    65  66  75  74;
    66  67  76  75;
    67  68  77  76;
    68  69  78  77;
    69  70  79  78;
    70  71  80  79;
    71  72  81  80;
];

E = 200000;
E1 = 1e-6;
v = 0.3;

L = 9;
A = 9;

nNodes = length(MeshNodes);
nEdges = 4;
nElems = length(MeshElems);
nNodeDofs = 2;

nElemDofs = nNodeDofs * 4;
nMeshDofs = nNodes * nNodeDofs;

D = E/(1 - v*v) * [1 v 0; v 1 0; 0 0 (1 - v)/2];
D1 = E1/(1 - v*v) * [1 v 0; v 1 0; 0 0 (1 - v)/2];

Forces = zeros(nMeshDofs, 1);

d = [1, 1 + nNodes, 2, 2 + nNodes, 3, 3 + nNodes, 4, 4 + nNodes, 5,...
    5 + nNodes, 6, 6 + nNodes, 7, 7 + nNodes, 8, 8 + nNodes, 9, ...
    9 + nNodes, 73+nNodes, 74+nNodes, 75+nNodes, 76+nNodes, 77+nNodes, ...
    78+nNodes, 79+nNodes, 80+nNodes, 81+nNodes];
dv = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1];

hold on;

function ShowBlue(MeshElems, MeshNodes, nElems)
  for i = 1:nElems
    if i ~= 28 && i~=29 && i~=30 && i ~= 20 && i~=21 && i~=22 && i~=27 && i~=42
      elem_nodes = [MeshElems(i, :), MeshElems(i, 1)];
      plot(MeshNodes(elem_nodes, 1), MeshNodes(elem_nodes, 2), 'b-o');
    end
  end
end

function ShowRed(MeshElems, MeshNodes, nElems)
  for i = 1:nElems
    if i ~= 28 && i~=29 && i~=30 && i ~= 20 && i~=21 && i~=22 && i~=27 && i~=42
      elem_nodes = [MeshElems(i, :), MeshElems(i, 1)];
      plot(MeshNodes(elem_nodes, 1), MeshNodes(elem_nodes, 2), 'r-o');
    end
  end
end

ShowBlue(MeshElems, MeshNodes, nElems);

function result = Nx(row, elemSize)
  y = row(2);
  a = elemSize(1);
  b = elemSize(2);
  result = 1/(a * b) * [-(b/2-y);
                        (b/2-y);
                        (b/2+y);
                        -(b/2+y)];
end

function result = Ny(row, elemSize)
  x = row(1);
  a = elemSize(1);
  b = elemSize(2);
  result = 1/(a * b) * [-(a/2-x);
                        -(a/2+x);
                        (a/2+x);
                        (a/2-x)];
end

function result = B(row, elemSize)
  nx = Nx(row, elemSize);
  ny = Ny(row, elemSize);
  result = [
    nx(1), nx(2), nx(3), nx(4), 0,     0,     0,     0;
    0,     0,     0,     0,     ny(1), ny(2), ny(3), ny(4);
    ny(1), ny(2), ny(3), ny(4), nx(1), nx(2), nx(3), nx(4);
  ];
end

function result = Ke(nElemDofs, D, elemSize)
  Ke = zeros(nElemDofs, nElemDofs);
  a = elemSize(1);
  b = elemSize(2);
  gp = [-1/sqrt(3) * a/2, 1/sqrt(3) * b/2];
  gw = [1 * a/2, 1 * b/2];
  for i = 1:2
    xi = gp(i);
    wi = gw(i);
    for j = 1:2
      eta = gp(j);
      wj = gw(j);
      Bm = B([xi, eta], elemSize);
      Ke = Ke + (Bm.' * D * Bm) * (wi * wj);
    end
  end
  result = Ke;
end


K_global = zeros(nMeshDofs, nMeshDofs);
for i = 1:size(MeshElems,1)
  elem = MeshElems(i,:);
  if i ~= 28 && i~=29 && i~=30 && i ~= 20 && i~=21 && i~=22 && i~=27 && i~=42
    Ke_local = Ke(nElemDofs, D, elemSize);
  else
    Ke_local = Ke(nElemDofs, D, elemSize);
  endif
  for k1 = 0:nNodeDofs-1
    for k2 = 0:nNodeDofs-1
      for row_loc = 1:nEdges
        row_glb = elem(row_loc);
        for col_loc = 1:nEdges
          col_glb = elem(col_loc);
          K_global(row_glb + k1*nNodes, col_glb + k2*nNodes) = ...
            K_global(row_glb + k1*nNodes, col_glb + k2*nNodes) + ...
            Ke_local(row_loc + k1*nEdges, col_loc + k2*nEdges);
        end
      end
    end
  end
end

function result = Pi(F, K, d, dv)
  for i = 1:length(F)
    for j = 1:length(d)
      k = d(j);
      if i ~= k
        F(i) = F(i) - K(i, k)*dv(j);
      end
    end
  end
  result = F;
end

function result = Pk(F, K, d, dv)
  for j = 1:length(d)
    k = d(j);
    F(k) = K(k, k)*dv(j);
  end
  result = F;
end

Forces = Pi(Forces, K_global, d, dv);
Forces = Pk(Forces, K_global, d, dv);

for i = 1:length(d)
  k = d(i);
  tmp = K_global(k, k);
  K_global(:, k) = 0;
  K_global(k, :) = 0;
  K_global(k, k) = tmp;
end

solu = linsolve(K_global, Forces);

function result = de(solu, elem, nNodes)
    result = [solu(elem); solu(elem + nNodes)]';
end

function sigma = sigmaElemTL(elemIndex, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
  a = elemSize(1);
  b = elemSize(2);
  elem = MeshElems(elemIndex,:);
  de = de(solu, elem, nNodes)
  xi = -a/2; eta = b/2;
  B = B([xi, eta], elemSize);
  eps_local = B * de';
  sigma = D * eps_local;
end

function sigma = sigmaElemTR(elemIndex, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
  a = elemSize(1);
  b = elemSize(2);
  elem = MeshElems(elemIndex,:);
  de = de(solu, elem, nNodes)
  xi = a/2; eta = b/2;
  B = B([xi, eta], elemSize);
  eps_local = B * de';
  sigma = D * eps_local;
end


sigma1l = sigmaElemTL(57, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma2l = sigmaElemTL(58, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma3l = sigmaElemTL(59, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma4l = sigmaElemTL(60, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma5l = sigmaElemTL(61, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma6l = sigmaElemTL(62, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma7l = sigmaElemTL(63, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma8l = sigmaElemTL(64, solu, MeshElems, MeshNodes, nNodes, elemSize, D)

sigma1r = sigmaElemTR(57, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma2r = sigmaElemTR(58, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma3r = sigmaElemTR(59, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma4r = sigmaElemTR(60, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma5r = sigmaElemTR(61, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma6r = sigmaElemTR(62, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma7r = sigmaElemTR(63, solu, MeshElems, MeshNodes, nNodes, elemSize, D)
sigma8r = sigmaElemTR(64, solu, MeshElems, MeshNodes, nNodes, elemSize, D)

ms1 = (sigma1l + sigma1r) / 2;
ms2 = (sigma2l + sigma2r) / 2;
ms3 = (sigma3l + sigma3r) / 2;
ms4 = (sigma4l + sigma4r) / 2;
ms5 = (sigma5l + sigma5r) / 2;
ms6 = (sigma6l + sigma6r) / 2;
ms7 = (sigma7l + sigma7r) / 2;
ms8 = (sigma8l + sigma8r) / 2;
ms_avg = ms1 + ms2 + ms3 + ms4 + ms5 + ms6 + ms7 + ms8

msM = [ ms_avg(1), ms_avg(3);
        ms_avg(3), ms_avg(2) ];

msV = (msM * [0;1]) * a;
Eyy = (msV(2) * 8 * b) / (8 * a * dv(end));

scale = 1;
NodesD = MeshNodes + scale * [solu(1:nNodes,:), solu(nNodes+1:2*nNodes,:)];

ShowRed(MeshElems, NodesD, nElems);
hold off;
