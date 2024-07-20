function [numInt, Int1, Int2] = Intersection_Of_Vector_With_Ellipsoid(P, v, a, b, c)
    % Intersection_Of_Vector_With_Ellipsoid computes the intersections, Int, of a
    % vector and a canonical ellipsoid with half-axes lengths a, b, c.

    syms t
    L = P+(t*v);
    x = L(1); % P(1)+(t*v(1)); 
    y = L(2); % P(2)+(t*v(2));
    z = L(3); % P(3)+(t*v(3));

    equ = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 == 0;
    t_solved = solve(equ);

    if isreal(t_solved) % Check if any intersections occur
        Int1 = double([subs(L(1), t, t_solved(1)), subs(L(2), t, t_solved(1)), subs(L(3), t, t_solved(1))]);
        Int2 = double([subs(L(1), t, t_solved(2)), subs(L(2), t, t_solved(2)), subs(L(3), t, t_solved(2))]);
        if Int1 == Int2
            numInt = 1;
        else 
            numInt = 2;
        end
    else % No intersections if t_solved contains imaginary numbers
        Int1 = [NaN,NaN,NaN];
        Int2 = [NaN,NaN,NaN];
        numInt = 0;
    end
end