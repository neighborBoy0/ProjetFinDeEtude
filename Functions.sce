m = mode();
mode(-1);


// the position error between the current and final positions for the end effector
function result = F1(CPosition, FPosition)
    try
        positions = [CPosition; FPosition];
        result = nan_pdist(positions, 'euclidean');  // calcule euclidean distance between CPosition and FPosition
    catch
        [error_message,error_number]=lasterror(%t)
        disp("There is an error in function F1, error message:" + error_message);
    end
endfunction

// The purpose of this term is to keep the control points (the center of the joints in this case) away from the various obstacles
//
// M is the number of obstacles
// N is the total number of control points
function result = F2(M, N)
    try
        f2 = 0
        // calcule sum of the euclidean distance between each joint and each obstacles 
        for i=1:N
            for j=1:M
                positions = [robot.links(i).r; [0.5, 0.5, 0.5]]; // Suppose there is only one obstacle at 0.5, 0.5, 0.5
                f2 = f2 + nan_pdist(positions, 'euclidean');
            end
        end
        result = f2;
    catch
        [error_message,error_number]=lasterror(%t)
        disp("There is an error in function F2, error message:" + error_message);
    end
endfunction

// The purpose of this term is to maximize the manipulability of the robot
//
// Q: 
function result = F3(Q)
    // Q = []
    // loop = size(robot.links);        // the number of joints
    // for i=1:loop
        // T = r2t(robot.links(6).I);
        // Qtemp = tr2q(T1)
        // Q = [Q, Qtemp(2)]
    // end
    J = jacob0(robot, Q, 'rot')         // it was jacobn
    result = sqrt(det(J * J'));
endfunction

function F4()
    
endfunction