% Author:  Bikash Kunwar/ bzk0067@auburn.edu
% Date: 2024-09-17
% Assignment Name: hw01

classdef hw01
    methods (Static)

        function p1()
            % This function only contains comments. Fill the following table. Do not write any code here.
            % :return: no returns

            % Write your result and explanation for each command here.
            %
            % commands         |  results      | explanations
            % -----------------|---------------|-----------------------------------
            % eps              |   2.2204e-16            |
            % realmax          | 1.7977e+308              |
            % realmin          | 2.2251e-308              |
            % 1 + eps - 1      |2.2204e-16               |
            % 1 + eps/2 - 1    |   0            |
            % realmin/1e10     |  2.2251e-318             |
            % realmin/1e16     |  0             |
            % realmax*10       |    Inf           |
        end

        function s_n = p2(n, choice)
            % This function computes the Archimedes' method for pi.
            % :param n: the number of sides of the polygon
            % :param choice: 1 or 2, the formula to use
            % :return: s_n, the approximation of pi using Archimedes' method.

            % Tabulate the error of |s_n - pi| for n = 0, 1, 2, ..., 15 and choice = 1 and 2.
            % for both choices of formulas.
            % n     | choice 1 | choice 2
            % ------|----------|----------
            %  0      0.322508961547962       0.322508961547962
            %  1        0.0737976555836819      0.0737976555836797
            %  2        0.0180672885077167      0.0180672885077082
            %  3       0.00449356154167369      0.0044935615416426
            %  4       0.00112194605578031     0.00112194605557647
            %  5      0.000280396390072823    0.000280396390031967
            %  6      7.00934652750895e-05    7.00934670572195e-05
            %  7      1.75230097285706e-05    1.75230148982131e-05
            %  8      4.38073354391832e-06    4.38073173514297e-06
            %  9       1.0952270628195e-06    1.09518156143906e-06
            % 10      2.74283840084877e-07    2.73795305538727e-07
            % 11      7.20327983927405e-08    6.84488226099234e-08
            % 12       1.8151752101403e-08    1.71122067627039e-08
            % 13      3.46889068580936e-08     4.2780530229436e-09
            % 14       1.8151752101403e-08    1.06951514311504e-09
            % 15      7.17707820285796e-07    2.67380784180205e-10

            % Explanation of the results (why there is a difference between the two choices):
            % As n increases, pn becomes incresingly closer to zero, and sqrt(1+pn^2)
            % tends to 1.
            % The absolute error when using formula 1 is greater than when using formula 2. This
            % is because of the cancellation error that arises when subtracting two
            % almost equal values: (sqrt(1+pn^2) and 1
            % Since formula 2 uses no such subtraaction, its error decreases
            % monotonical and then reamins constant.

            % Write your code here
            % n = 0,1,2,....
            % ind = 1,2,3,...
            p0 = 1/sqrt(3);
            p(1) = p0; % n = 0
            if choice == 1
                % Use the 1st formula
                % s_n = inf; % Write your code here




                %% formula 1

                for n_counter = 0:n
                    ind = n_counter+1;
                    sqrt_term(ind) = sqrt(1+p(ind)^2);
                    p(ind+1) = (sqrt_term(ind)-1)/p(ind);
                    S(ind) = 2^n_counter*6*p(ind);
                end
                s_n = S(end);
            else
                % Use the 2nd formula
                % s_n = inf; % Write your code here
                for n_counter = 0:n
                    ind = n_counter+1;
                    sqrt_term(ind) = sqrt(1+p(ind)^2);
                    p(ind+1) = p(ind) /(sqrt_term(ind) +1);
                    S(ind) = 2^n_counter*6*p(ind);
                end
                s_n = S(end);
            end

        end

        function s = p3(a)
            % This function computes the Kahan summation algorithm.
            % :param a: a vector of numbers
            % :return: summation of the vector a using Kahan summation algorithm

            % s = inf; % Write your code here

            n = length(a);

            %% Initialization
            j = 1;
            e(1) = 0;
            S(1) = a(1);

            %% Calculation
            while j<n
                j = j+1;
                y(j) = a(j)-e(j-1);
                S(j) = S(j-1) + y(j);
                e(j) = (S(j)-S(j-1)) - y(j);
            end
            s = S(end);
        end

        function p4(a)
            % This function test the performance of Kahan summation algorithm against native sum.
            % :param a: a vector of numbers in double precision.
            % :return: no returns

            % Test this function with a = rand(n, 1) with various size n multiple times.
            % Summarize your findings below.
            %
            % Findings:
            %
            % Among different trial set of numbers in x, the Kahan summation and the
            % built in summation function sum() under single precision gave the same
            % answer for most cases and a slightly different for ver few
            % cases. Cound not tell which one is more accurate, although
            % Khar summation should be given matlab's sum suffers floating
            % point roundoff errors. 

            single_a = single(a); % Convert a to single precision
            s = hw01.p3(a); % Kahan sum of a under double precision (regarded as truth).

            single_Kahan_s = hw01.p3(single_a); % Kahan sum of single_a under single precision.
            single_naive_s = sum(single_a); % Naive sum of single_a under single precision.

            disp(['Error of naive sum under single precision: ', num2str(single_naive_s-s)]);
            disp(['Error of Kahan sum under single precision: ', num2str(single_Kahan_s-s)]);
        end

        function s = p5(a)
            % For 6630.
            % This function computes summation of a vector using pairwise summation.
            % :param a: a vector of numbers
            % :return: the summation of the vector a using pairwise summation algorithm.

            % ! You may need to create a helper function if your code uses recursion.

            % Rewrite the p4 test function to test this summation method.
            % Summarize your findings below.
            %
            % Findings:
            %
            %
            %
            %
            %


            s = inf; % Write your code here.

        end
    end
end