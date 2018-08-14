function folder = datafolder()

    if strfind(getenv('OS'),'Windows')
        folder = fullfile('C:/Users/Dan/Documents/Box Sync/COHCON_DATA');
    else
        folder = fullfile('~/Box Sync/COHCON_DATA');
    end
