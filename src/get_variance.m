function state_var = get_variance(P,N_col_states,N_states)

    var_py = diag(P(1:N_col_states:N_states,1:N_col_states:N_states)) ...
        + diag(P(3:N_col_states:N_states,3:N_col_states:N_states)) ...
        + diag(P(5:N_col_states:N_states,5:N_col_states:N_states)) ...
        + diag(P(7:N_col_states:N_states,7:N_col_states:N_states)) ...
        + diag(P(11:N_col_states:N_states,11:N_col_states:N_states)) ...
        + 2*diag(P(1:N_col_states:N_states,3:N_col_states:N_states)) ...
        + 2*diag(P(1:N_col_states:N_states,5:N_col_states:N_states)) ...
        + 2*diag(P(1:N_col_states:N_states,7:N_col_states:N_states)) ...
        + 2*diag(P(1:N_col_states:N_states,11:N_col_states:N_states)) ...
        + 2*diag(P(3:N_col_states:N_states,5:N_col_states:N_states)) ...
        + 2*diag(P(3:N_col_states:N_states,7:N_col_states:N_states)) ...
        + 2*diag(P(3:N_col_states:N_states,11:N_col_states:N_states)) ...
        + 2*diag(P(5:N_col_states:N_states,7:N_col_states:N_states)) ...
        + 2*diag(P(5:N_col_states:N_states,11:N_col_states:N_states)) ...
        + 2*diag(P(7:N_col_states:N_states,11:N_col_states:N_states));
    
    var_i = diag(P(9:N_col_states:N_states,9:N_col_states:N_states));
    var_e = diag(P(13:N_col_states:N_states,13:N_col_states:N_states));
    
    state_var = zeros(12,1);
    state_var(1:3:end) = var_py;
    state_var(2:3:end) = var_e;
    state_var(3:3:end) = var_i;


end

