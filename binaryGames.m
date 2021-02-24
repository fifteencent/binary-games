keywords = ["merlin", "lights out"];
keywords = whichProgramToRun(keywords);
if keywords(1,2) == "lights out"
    prompt = 'Enter your Lights Out binary matrix (1 is on, 0 is off) in the format [r1v1, r1v2; r2v1, etc.]: ';
    Ai = input(prompt);
    Ai = isInputCorrect(Ai, prompt)
    [dim,~] = size(Ai);
    prompt = 'Enter your final Lights Out pattern: ';
    Af = input(prompt);
    Af = isInputCorrect(Af, prompt)
    soln = lightsOut(Ai(:), Af(:),dim)
    moves = sum(soln,'all');
    disp("with " + string(moves) + " moves")
end

function soln = lightsOut(Ai, Af, dim)
    B = mod(Af - Ai, 2);
    b = B(:);
    M = createM(dim);
    
    [~,f] = size(M);
    [~,h] = size(b);

    C = [M,b];
    C = rowEchelon(C);
    D = C(:,1:f);
    B = C(:,f + 1:f + h);
    rank = findRankOfEchelon(D);
    if rank < findRankOfEchelon(B)
        disp('There is no solution.');
    else 
        nD = D(1:rank,1:rank);
        B = B(1:rank,1:min([rank,h]));
        soln = mod(nD\B,2);
        null = nullRows(D);
        [a,~] = size(null);
        newD = nullMatrix(nD,null);
        newnull = zeros(f,a);
        for i = 1:a
            newB = mod(-1*D(:,null(i)),2);
            newnull(:,i) = recreate(mod(newD\newB(1:rank), 2),null,null(i));
        end
        delRows = zeros(dim^2-rank,1);
        soln = [soln;delRows];
        soln = smallestSoln(soln, newnull,dim);
        soln = reshape(soln, dim, dim);
    end
end

function new = recreate(A,null,index)
    [m,~] = size(null);
    [n,~] = size(A);
    new = zeros(m+n,1);
    for i = 1:m+n
        if ismember(i, null)
            new(i) = 0 + (i == index);
        else
            new(i) = A(1);
            A = A(2:end);
        end
    end
end

function nullD = nullMatrix(A,others)
    nullD = [];
    [m,~] = size(others);
    others = [0; others];
    for n=2:m
        nullD = [nullD, A(:,others(n-1)+1:others(n)-1)];
    end
    nullD = [nullD A(:,others(m)+1:end)];
end

%TODO: fix indexing issue
function soln = smallestSoln(A, nullM,dim)
    A = mod(A,2);
    [~,numNulls] = size(nullM);
    if isempty(nullM)
        soln = A;
    else
        nullM = [zeros(dim^2,1),nullM];
        solns = nullM;
        solns(:,1) = A;
        counts = zeros(numNulls+1,1);
        counts(1) = sum(A,'all');
        for i = 2:numNulls+1
            null1 = nullM(:,i);
            newTemp = smallestSoln(A + null1,nullM(:,i+1:end),dim);
            solns(:,i) = newTemp;
            counts(i) = sum(newTemp,'all');
        end
        cmin = min(counts);
        imin = find(counts == cmin);
        [a,~] = size(imin);
        if a > 1
            imin = imin(1);
        end
        soln = mod(solns(:,imin),2);
    end
end

function M = createM(dim)
    M = zeros(dim^2);
    i = 1;
    for c = 1:dim
        for r = 1:dim
            M(:,i) = createMCols(r,c,dim);
            i = i + 1;
        end
    end
end

function mcol = createMCols(row,col, dim)
    mcol = zeros(dim);
    minR = max([row-1,1]);
    maxR = min([row+1,dim]);
    for r = minR:maxR
        mcol(r,col) = 1;
    end
    minC = max([col-1,1]);
    maxC = min([col+1,dim]);
    mcol(row,minC) = 1;
    mcol(row,maxC) = 1;
    mcol = mcol(:);
end

function null = nullRows(A)
    [dim,~] = size(A);
    rowCounts = zeros(dim,1);
    rowCounts(1,1) = 1;
    null = zeros(dim, 1);
    for n = 2:dim
        rowCounts(n) = findFirstNonZeroRow(transpose(A(n,:)), 1);
        if rowCounts(n) < 0
            rowCounts(n) = dim;
        end
        count = rowCounts(n);
        prevCount = rowCounts(n-1);
        diff = count - prevCount;
        if diff > 1
            for i = 1:diff-1
                null(prevCount+i) = prevCount+i;
            end
        elseif diff == 0
             null(count) = count;
        end
    end
    null = nonzeros(null);
end

function x = rowEchelon(A)
    [m,n] = size(A);
    if m > n
        m = n;
    end
    
    i = 1;
    while i <= m
        A = mod(eliminateColumn(A,i),2);
        i = i + 1;
    end
    x = A;
end

%A is the matrix, d is the column
function x = eliminateColumn(A,d)
    a = findFirstNonZeroRow(A, d);
    [m,~] = size(A);
    if a > d
        b = A(d,:);
        A(d,:) = A(a,:);
        A(a,:) = b;
    end
    while m > d
        if A(m,d) ~= 0 && A(d,d) ~= 0
            c = A(m,d)/A(d,d);
            A(m,:) = mod(A(m,:) - c * A(d,:),2);
        end
        m = m - 1;
    end
    x = A;
end

function x = findFirstNonZeroRow(A, d)
    [m,~] = size(A);
    i = d;
    % keeps counting until the index does not equal 0
    while A(i,d) == 0
        i = i + 1;
        if i > m
            break
        end
    end
    if i > m
        i = -1;
    end
    x = i;
end

function x = findRankOfEchelon(A)
    [m,n] = size(A);
    if m > n
        m = n;
    end
    
    %keeps counting until the index is not equal to 0
    while sum(A(m,:)) == 0 && m > 0
        m = m - 1;
        if m < 1
            break
        end
    end
    x = m;
end

function x = isInputCorrect(y, prompt)
    prompt2 = 'Is this what you entered (this is a yes or no question)? ';
    disp(y)
    str = input(prompt2, 's');
    if contains(str, 'n')
        x = input(prompt);
        isInputCorrect(x, prompt);
    elseif contains(str,'y')
        x = y;
    else
        disp('This is a yes or no question! (Please answer yes or now)');
        x = isInputCorrect(y,prompt);
    end
end

function keys = whichProgramToRun(keywords)
    prompt = 'What would you like to find? ';
    str = input(prompt, 's');
    keys = strings(1,length(keywords));
    options = "";
    for n = 1:length(keywords)
        word = keywords(n);
        if contains(word, str)
            keys(n) = word;
        else
            keys(n) = ' ';
        end
        options = options + '"' + word + '",  ';
    end
    options = char(options);
    options = options(1:end-3);
    if ~max(contains(keywords, str))
        disp("Please use one of the keywords: " + options + "!")
        keys = whichProgramToRun(keywords);
    end
end