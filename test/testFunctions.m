function testFunctions(n)
    % the result shows that reading from obj is 10 times slower than from
    % normal variable, and writing is even slower.
    %
    % The differences between x and opt.x are negligible
    z=zeros(100,1);
    tic;
    x=z;
    for i=1:n
        x=x+1;
    end
    scriptTime=toc;

    tic;
    add(n);
    innerfuctionTime=toc;

    tic;
    addz(z,n);
    functoinTime=toc;
    display([scriptTime, innerfuctionTime, functoinTime]);


    test1=@test;

    a=10;
    b=20;
    display([a b]);

    test(30)

    display([a b]);

    test2=@test;

    test1(40)

    test2(50)

    display([a b]);

    function add(n)
        x=z;
        for i=1:n
            x=x+1;
        end
    end
    function test(b)
        display([a b]);
        a=a-5;
        b=b-5;
        display([a b]);
    end
end

function addz(z,n)
    x=z;
    for i=1:n
        x=x+1;
    end
end
