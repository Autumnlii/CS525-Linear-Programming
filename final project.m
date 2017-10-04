
%% Problem 1.
%Name: Qiuying Li
%Time: 4/28/2016
% Write code in Matlab and CVX to solve the QP formulation above.
%%
% (a) Test your routine on the training data obtained by setting
% fracTest=0.1 and reord=0 in help wdbcData, using the value ?= 0.001.
% (This will yield a training data set with 512 data points, consisting of
% records 1 through 512 from the file wdbc.data.)

clear;clc;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.1,0);

%Set up two matrix, M and B.
%M has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
B=zeros(ntrain,30);
M=zeros(ntrain,30);
k=0;
m=0;
% Then split training set into the two matrices M and B
i=1;
while (i <= ntrain)
    if train(i,1) ~= 1
         k = k+1;
        B(k,:) = train(i,2:31);
    else
        m = m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
B = B(1:k,:);
M = M(1:m,:);
% Set up the CVX configuration
% We have four variables, Omega which 
% omega is the normal vector of the hyperplane
% Gamma is the offset of the hyperplane away from the origin
mu = 0.001;
em = ones(m,1);
ek = ones(k,1);

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k)
minimize( (em'*y)/m + (ek'*z)/k +  (mu/2)*(Omega'*Omega))  
subject to
M*Omega - em*Gamma + y >= em
-B*Omega + ek*Gamma + z >= ek
z >= 0
y >= 0
cvx_end

% Test the model on training set
countmisplaced = 0;
test = zeros(ntrain,1);
j=1;
while (j <= ntrain)
    test(j) = train(j,2:31)*Omega - Gamma;
     if ( test(j) < 0 && train(j,1)==1 )
        countmisplaced = countmisplaced + 1;
     end
     if ( test(j) > 0 && train(j,1)==0 )
        countmisplaced = countmisplaced + 1;
     end   
    j=j+1;
end
%Display all the results
disp('Omega is:');
disp(Omega);
fprintf('Offset from the origin   Gamma = %d\n',Gamma);
fprintf('Number of misclassified points = %d out of %d\n',countmisplaced,ntrain);
fprintf('Optimal objective = %d\n',cvx_optval);

%%Problem 1
% (b) Test your routine on the training data obtained by setting
% fracTest=0.15 and reord=1 in help wdbcData, using ?= 0.001.
% (This will yield a training matrix of 484 records randomly selected from
% the 569 samples in the file wdbc.data.)
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.15,1);
%Similarly, set up two matrices, M and B
%M has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
M=zeros(ntrain,30);
B=zeros(ntrain,30);
m = 0;
k = 0;
% Split training set into the two matrices M and B
i=1;
while (i<=ntrain)
    if train(i,1) ~= 1
         k = k+1;
        B(k,:) = train(i,2:31);
    else
       m = m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
B = B(1:k,:);
M = M(1:m,:);
%There are four variables Omega(normal vector of hyperplane) , Gamma(offset of hyperplane), y and z
em = ones(m,1);
ek = ones(k,1);
mu = 0.001;

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k),
minimize((em'*y)/m+(ek'*z)/k+(mu/2)*(Omega'*Omega))
subject to
M*Omega - em*Gamma + y >= em
-B*Omega + ek*Gamma + z >= ek
y >= 0
z >= 0
cvx_end

%Test training set
countmisplaced = 0;
test = zeros(ntrain,1);
j=1;
while(j<=ntrain)
    test(j) = train(j,2:31)*Omega-Gamma;
    if (test(j)<0&&train(j,1)==1)
        countmisplaced=countmisplaced + 1;
    end
    if (test(j)>0&&train(j,1)==0)
        countmisplaced=countmisplaced + 1;
    end
    j=j+1;
end
%Display the results
disp('Omega is:');
disp(Omega);
fprintf('Offset from the origin   Gamma = %d\n',Gamma);
fprintf('Number of misclassified points = %d out of %d\n',countmisplaced,ntrain);
fprintf('Optimal objective = %d\n',cvx_optval);

%% 2.
% Modify your code for part 1, write a program to obtain the separating
%plane on the training set, and determine the number of misclassified points on the corresponding
% testing set, for the following cases (use ?= 0.001 for each):
% (a) fracTest=0.1 and reord=0
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.1,0);
% Similarly, set training set into M and B.
% M has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
M=zeros(ntrain,30);
B=zeros(ntrain,30);
m =0;
k =0;
i=1;
while (i<=ntrain)
    if train(i,1) ~= 1
        k=k+1;
        B(k,:) = train(i,2:31);
    else
        m=m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
M = M(1:m,:);
B = B(1:k,:);
%CVX configuration
%four variables Omega (same as before), Gamma,, y, z
mu = 0.001;
em = ones(m,1);
ek = ones(k,1);

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k)
minimize((em'*y)/m+(ek'*z)/k+(mu/2)*(Omega'*Omega))
subject to
M*Omega -em*Gamma + y >= em
-B*Omega +ek*Gamma + z >= ek
y>= 0
z>= 0
cvx_end
% Test
countmisplaced= 0;
j=1;
while (j<=ntest)
    test_t(j) = test(j,2:31)*Omega - Gamma;
    if (test_t(j)>0 && test(j,1)==0 )
        countmisplaced = countmisplaced+1;
    end
    if (test_t(j)<0 && test(j,1)==1 )
        countmisplaced = countmisplaced+1;
    end
    j=j+1;
end
fprintf('When fracTest = 0.1 and reord = 0\n')
fprintf('There are %d misclassified points.\n',countmisplaced)

% (b) fracTest=0.15 and reord=0
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.15,0);
%Split into M and B.
%M has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
M=zeros(ntrain,30);
B=zeros(ntrain,30);
m =0;
k =0;
i=1;
while (i <=ntrain)
    if train(i,1) ~= 1
         k = k+1;
        B(k,:) = train(i,2:31);
    else
       m = m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
M = M(1:m,:);
B = B(1:k,:);
%CVX configuration
%four variables, Omega, Gamma, y and z (same as before)
em = ones(m,1);
ek = ones(k,1);
mu = 0.001;

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k)
minimize((em'*y)/m+(ek'*z)/k+(mu/2)*(Omega'*Omega))
subject to
M*Omega - em*Gamma + y >= em
-B*Omega + ek*Gamma + z >= ek
y>= 0
z>= 0
cvx_end
% Test
countmisplaced = 0;
test_t= zeros(ntest,1);
j=1;
while (j<=ntest)
    test_t(j) = test(j,2:31)*Omega - Gamma;
    if (test_t(j)>0 && test(j,1)==0 )
        countmisplaced = countmisplaced+1;
    end
    if (test_t(j)<0 && test(j,1)==1 )
        countmisplaced = countmisplaced+1;
    end
    j=j+1;
end
fprintf('When fracTest = 0.15 and reord = 0\n')
fprintf('There are %d misclassified points.\n',countmisplaced)

% (c) fracTest=0.20 and reord=1
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.2,1);

%Split into M and B where
%M has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
M=zeros(ntrain,30);
B=zeros(ntrain,30);
m = 0;
k = 0;
i=1;
while (i <=ntrain)
    if train(i,1) ~= 1
         k = k+1;
        B(k,:) = train(i,2:31);
    else
       m = m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
B = B(1:k,:);
M = M(1:m,:);
% CVX configuration
%we have four variables, Omega, Gamma, y and z (same as before)
mu = 0.001;
em = ones(m,1);
ek = ones(k,1);

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k)
minimize((em'*y)/m +(ek'*z)/k +(mu/2)*(Omega'*Omega))
subject to
M*Omega - em*Gamma + y >= em
-B*Omega+ ek*Gamma + z >= ek
y>= 0
z>= 0
cvx_end

%test
countmisplaced = 0;
test_t = zeros(ntest,1);
j=1;
while (j <=ntest)
    test_t(j) = test(j,2:31)*Omega - Gamma;
    if (test_t(j)<0 && test(j,1)==1 )
        countmisplaced=countmisplaced + 1;
    end
    if (test_t(j)>0 && test(j,1)==0 )
        countmisplaced=countmisplaced + 1;
    end
    j=j+1;
end
fprintf('When fracTest is 0.20 and reord = 1\n')
fprintf('There are %d misclassified points.\n',countmisplaced)

% (d)when fracTest=0.05 and reord=1
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.05,1);

%Similarly, split training set into M and B whereM has all information for malignant cells in the training set, while B
%has information for benign cells in the training set.
B=zeros(ntrain,30);
M=zeros(ntrain,30);
m= 0;
k= 0;
i=1;
while (i <=ntrain)
    if train(i,1) ~= 1
         k = k+1;
        B(k,:) = train(i,2:31);
    else
       m = m+1;
        M(m,:) = train(i,2:31);
    end
    i=i+1;
end
M = M(1:m,:);
B = B(1:k,:);
%CVX configuration
%Omega, Gamma, y and z are the four variables (same as before)
em = ones(m,1);
mu = 0.001;
ek = ones(k,1);

cvx_begin quiet
variables Omega(30) Gamma(1) y(m) z(k)
minimize((em'*y)/m+(ek'*z)/k +(mu/2)*(Omega'*Omega))
subject to
M*Omega - em*Gamma + y >= em
-B*Omega + ek*Gamma + z >= ek
y>= 0
z>= 0
cvx_end

%test
countmisplaced = 0;
test_t = zeros(ntest,1);
j=1;
while (j <=ntest)
    test_t(j) = test(j,2:31)*Omega - Gamma;
    if (test_t(j) < 0 && test(j,1)==1 )
        countmisplaced=countmisplaced + 1;
    end
    if (test_t(j)>0 && test(j,1)==0 )
        countmisplaced=countmisplaced + 1;
    end
    j=j+1;
end
fprintf('When fracTest is 0.05 and reord = 1\n')
fprintf('there are %d misclassified points.\n',countmisplaced)

%%Problem 3
clear;
[train,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.12,0);
%reset the count for m and b.
countBtype=0;
countMtype=0;
j=1;
while (j<= ntrain)
    %if not equal to 1, it is b type
    if train(j,1) ~= 1
        countBtype = countBtype+1;
    else
        countMtype = countMtype+1;
    end
    j=j+1;
end
%Put result in array to store
array_B = zeros(countBtype,2);
array_M = zeros(countMtype,2);
Count = ntrain;
m=2;
while (m <=30)
    firstfeature= train(:,m);
    i=m+1;
    while (i <=31)
        secondfeature =train(:,i);
        l1= 0;
        l2= 0;
        k=1;
        while (k<=ntrain)
            if train(k,1) ~= 1
                l2 = l2+1;
                array_B(l2,1) = firstfeature(k,1);
                array_B(l2,2) = secondfeature(k,1);
            else
                l1 = l1+1;
                array_M(l1,1) = firstfeature(k,1);
                array_M(l1,2) = secondfeature(k,1);
            end
            k=k+1;
        end
        %CVX configuration when miu=0.0008
        %Four variables Omega, Gamma, y and z
        em = ones(countMtype,1);
        ek = ones(countBtype,1);
        miu = 0.0008;
        %begin CVX
        cvx_begin quiet
        variables Omega(2) Gamma(1) y(countMtype) z(countBtype)
        minimize((em'*y)/countMtype+(ek'*z)/countBtype+(miu/2)*(Omega'*Omega))
        subject to
        array_M*Omega - em*Gamma +y >= em
        -array_B*Omega + ek*Gamma +z >= ek
        y>= 0
        z>= 0
        cvx_end
        % Testing
        countmisplaced = 0;
        test_t = zeros(ntrain,1);
        j=1;
        while (j <=ntrain)
            test_t(j) = train(j,[m i])*Omega - Gamma;
            if (test_t(j)>0 && train(j,1)==0 )
                countmisplaced = countmisplaced + 1;
            end
            if (test_t(j)<0 && train(j,1)==1 )
                countmisplaced=countmisplaced + 1;
            end
            j=j+1;
        end
        %discuss and display the results
        if countmisplaced < Count
            fprintf('attributes %2d %2d: misclass %3d\n',m,i, countmisplaced);
            Count=countmisplaced;
            best_omega = Omega;
        end
        i=i+1;
    end
    m=m+2;
end
% Problem 4
%After the conclusion of problem 3, we have known that the best case is when m = 15, i = 29
%testing
test_t = zeros(ntest,1);
countmisplaced = 0;
j=1;
while (j <=ntest)
    test_t(j) = test(j,[15 29])*best_omega - Gamma;
    if ( test_t(j) > 0 && test(j,1)==0 )
        countmisplaced=countmisplaced + 1;
    end
    if (test_t(j)<0 && test(j,1)==1)
        countmisplaced = countmisplaced + 1;
    end
    j=j+1;
end
fprintf('number of misclassified elements is');
disp(countmisplaced)