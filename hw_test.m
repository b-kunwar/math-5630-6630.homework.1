hw01_worker = hw01();

hw_assert(abs(hw01_worker.p2(0, 1) - 2 * sqrt(3)) < 1e-8);
hw_assert(abs(hw01_worker.p2(1, 1) - 12/sqrt(3)/(1 + sqrt(4/3)))<1e-8);

hw_assert(abs(hw01_worker.p2(0, 2) - 2 * sqrt(3))<1e-8);
hw_assert(abs(hw01_worker.p2(1, 2) - 12/sqrt(3)/(1 + sqrt(4/3)))<1e-8);

hw_assert(hw01_worker.p3([1, 2, 3]) == 6);
hw_assert(hw01_worker.p3([eps, -eps]) == 0);


n = 4;
hw01_worker.p4(rand(n,1));

function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end
