function findParameters(sino)
    [pw,np]=size(sino);
    sino=sino>=0.15;
    g1=sum(sino(1:pw/2,:))+pw/2;
    g2=sum(sino(pw/2+1:end,:))+pw/2;
    G1=fft(g1);
    G2=fft(g2);

    D = angle(G1./G2);


    
end
