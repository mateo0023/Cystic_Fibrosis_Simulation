from multiprocessing.pool import Pool
import pandas as pd
import numpy as np
import datetime
import zipfile
import sys
import os
import simulation
AirwayModel = simulation.AirwayModel

# We need to stop the simulations after any error, else data will be harder to analyze.
prev_np_err = np.seterr(all='raise')
# Will be the folder where all of the temporary files will be saved.
SUB_DIR = 'data'
# Will be the MAX files that can be on the SUB_DIR before they are sent to a ZIP file. Doesn't include extra-info
MAX_FILES = 5
# Max number of parallel processes.
MAX_PROCESSES = 8


# Run the simulation with the parameters in init_d and save its data. Will fill the complete pd.DataFrame
def runAll(init_d=None, save_extra=True, sub_dir=SUB_DIR):
    """
    Will create a AirwayModel object and fill up all of it. Once the model is complete it
    will save the results into a csv file. If an error is encountered the run will exit at that point.
    Optionally it will save extra info into a secondary CSV (if there are no errors):
        * Initial values of the variables.
        * Average values of the variables.
        * Runtime.
        * How many seconds each step represented.
        * Total steps the simulation took.
    If there is any error (independently of save_extra), the secondary CSV will contain the following:
        * Initial values of the variables.
        * Values for the step previous to the error.
        * The values of the variables when the error was raised.
        * Runtime.
        * How many seconds each step represented.
        * Total steps the simulation took.
        * And an error code with the step at which it was received.

    :param sub_dir: The directory to which to save the CSV files, defaults to SUB_DIR (global in the simulation).
    :param init_d: The initial values of the DataFrame. If None, it will still create the model, which will ask for the
            values individually.
    :param save_extra: Boolean, whether to save the CSV with the extra info. If there are any issues encountered,
            this will be overwritten to True.
    :return: The name of the file to which the data was saved (NOT the path, just the file name).
    """

    date_time = datetime.datetime.now()
    time = datetime.time(hour=date_time.hour, minute=date_time.minute, second=date_time.second)
    print('Starting run.\n\tTime:', time)
    # Create d_extr to save all the averages and runtime
    d_extr = {}
    if type(init_d) is AirwayModel:
        gen_data = init_d
    else:
        gen_data = AirwayModel(init_d)

    init_time = datetime.datetime.now()

    # Run all, but if an error is encountered, exit the run. It'll be easier to see the error in that case.
    # And save time if the error was caused only in one of the values due to initial conditions.
    if type(init_d) is AirwayModel:
        for step in range(1, len(gen_data)):
            # noinspection PyBroadException
            try:
                gen_data.fn_run(step)
            except:
                print('Unexpected error at step', step, '\n\t', sys.exc_info()[0:2])
                d_extr['err'] = [sys.exc_info()[0], step, None]
                # Drop the unfilled rows, they consume space and just complicate things.
                gen_data.drop(labels=range(step + 1, gen_data.max_steps))
                break
    else:
        for step in range(1, init_d['max_steps']):
            try:
                gen_data.fn_run(step)
            except:
                print('Unexpected error at step', step, '\n\t', sys.exc_info()[0:2])
                d_extr['err'] = [sys.exc_info()[0], step, None]
                # Drop the unfilled rows, they consume space and just complicate things.
                gen_data.drop(labels=range(step + 1, gen_data.max_steps))
                init_d['max_steps'] = step+1
                break
        # Turn all the activities back into concentration.
        for step in range(init_d['max_steps']):
            gen_data['aNa'][step] = gen_data.to_cons(gen_data['aNa'][step])
            gen_data['aCl'][step] = gen_data.to_cons(gen_data['aCl'][step])
            gen_data['aK'][step] = gen_data.to_cons(gen_data['aK'][step])

            gen_data['cNa'][step] = gen_data.to_cons(gen_data['cNa'][step])
            gen_data['cCl'][step] = gen_data.to_cons(gen_data['cCl'][step])
            gen_data['cK'][step] = gen_data.to_cons(gen_data['cK'][step])

    end_time = datetime.datetime.now()
    runtime = end_time - init_time

    file_name = 'Data_' + str(datetime.datetime.now())[:-5].replace(':', '-').replace(' ', '_')

    if gen_data.isCF:
        file_name += '_CF'
    else:
        file_name += '_NL'

    # Save the file
    if os.path.isdir(sub_dir):     # Checking if dir exists
        gen_data.to_csv(sub_dir + '/' + file_name + '.csv')
    else:
        gen_data.to_csv(file_name + '.csv')

    if 'err' in d_extr:
        # The step at which the error occurred
        step_err = d_extr['err'][1]
        # If exited from an error, save the relevant data: initial, final, and final-1
        for variable in (AirwayModel.variables[1::])[::-1]:
            d_extr[variable] = [gen_data[variable][0], gen_data[variable][step_err - 1],
                                gen_data[variable][step_err]]
        d_extr['runtime'] = [str(end_time - init_time), str(init_time), str(end_time)]
        d_extr['time_frame'] = [gen_data.time_frame, None, None]
        d_extr['max_step'] = [gen_data.max_steps, None, None]
        # Save the extra data about the error
        pd.DataFrame(d_extr, index=['Initial', 'Previous', 'Final (err)']).to_csv(
            sub_dir + '/' + file_name + '__error-info.csv')
    elif save_extra:
        # To ignore 'Time (min)' we start from the index 1
        for variable in AirwayModel.variables[1::]:
            av = np.sum(gen_data[variable]) / gen_data.max_steps
            d_extr[variable] = [gen_data[variable][0], gen_data[variable][round(gen_data.max_steps/2)],
                                gen_data[variable][gen_data.max_steps-1], av]
        d_extr['runtime'] = [str(end_time - init_time), None, None, None]
        d_extr['start-end'] = [str(init_time), str(end_time), None, None]
        d_extr['time_frame - max_step'] = [gen_data.time_frame, gen_data.max_steps, None, None]
        # Save the extra data about the run
        pd.DataFrame(d_extr, index=['initial', 'middle', 'final', 'average']).to_csv(
            sub_dir + '/' + file_name + '__extra-info.csv')

    date_time = datetime.datetime.now()
    time = datetime.time(hour=date_time.hour, minute=date_time.minute, second=date_time.second)
    del gen_data
    print('Done with the run.', '\n\tRuntime:', runtime, '\n\tAt time:', time)
    # noinspection PyUnboundLocalVariable
    return step


# Will return all the files inside a directory.
def flsInDir(directory='.'):
    """
    Will look at all the files in a directory and collect them into a list. Will walk using the OS Library
    :param directory: The path to the directory. Defaults to the current directory
    :return: A list with the paths to the files inside the directory. None if there are no files inside.
    """
    files_in = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            files_in.append(os.path.join(root, file))
    if len(files_in) == 0:
        return None
    else:
        return files_in


# Will delete the given folder, all of its contents and sub-contents.
def deleteFolder(directory=SUB_DIR):
    """
    Will delete the folder and all of it's contents, including sub-directories.
    :param directory: The directory to which to delete.
    :return: True if successful.
    """
    for root, dirs, files in os.walk(directory):
        for d in dirs:
            deleteFolder(d)
        for f in files:
            os.remove(os.path.join(root, f))
    os.rmdir(directory)
    return True


# Will save the file_name.csv and file_name_extra-info.csv files into one ZIP file.
def zipFiles(file_names=flsInDir(SUB_DIR), zip_name=None, delete_after=False):
    """
    Will create a zip file in the current directory which contains the files in the file_names
    :param file_names: Should be an iterable with the path and file names of the files to compress.
            Each item must be the path and name of the file.
            The files should have an extension, if a dot character ('.') is not detected in the 4th or 3rd position
            from the back, it will attach the file extension '.csv'
    :param zip_name: (string) The name of the ZIP file, it will add a '.zip' extension if it's not already in the name.
            It defaults to None, in which case it will generate one with 'Data_' and the current date.
    :param delete_after: boolean whether or not to delete the files after compressing. Defaults to False
    :return: True if completed successfully
    """
    if zip_name is None:
        zip_name = 'Data_' + str(datetime.datetime.now())[0:-10].replace(':', '-').replace(' ', '_') + '.zip'
    elif '.zip' not in zip_name:
        zip_name += '.zip'
    zip_file = zipfile.ZipFile(zip_name, 'w')

    for file in file_names:
        if '.' not in file[-4:-2]:
            try:
                zip_file.write(file + '.csv')
            except FileNotFoundError:
                print('Could not find the file', file + '.csv')
        else:
            try:
                zip_file.write(file)
            except FileNotFoundError:
                print('Could not find the file', file)

    print(zip_name, 'was created.')

    if delete_after:
        for file in file_names:
            try:
                # Check to see if there's a file format attached
                if '.' not in file[-4::-2]:
                    os.remove(file + '.csv')
                else:
                    os.remove(file)
            except FileNotFoundError:
                pass
        return True
    else:
        return True

def loadInitialValues():
    initial_values = []
    initial_pd = pd.read_csv('initial-values.csv')
    for i in range(len(initial_pd)):
        initial_values.append({})
        for var in initial_pd:
            if var == 'max_steps':
                initial_values[-1][var] = int(initial_pd[var][i])
            elif var == 'CF':
                initial_values[-1][var] = bool(initial_pd[var][i])
            else:
                initial_values[-1][var] = float(initial_pd[var][i])
    return initial_values

if __name__ == '__main__':
    # Empty the SUB_DIR folder before starting. If there are files inside, try to zip and delete them.
    if os.path.isdir(SUB_DIR):
        fls = flsInDir(SUB_DIR)
        if fls is not None:
            zipFiles(fls, 'OldData' + str(datetime.date.today()), delete_after=True)
        os.rmdir(SUB_DIR)

    # Try to load the 'initial-values.csv' file and run all of the simulations in parallel in sets of
    # MAX_PROCESSES using the multiprocess.Pool module.
    # If it can't be done (FileNotFound error), then run only one simulation, asking for all of the input values.
    try:
        # Prepare the initial values into a list of dictionaries.
        initial_values = loadInitialValues()

        print('Loaded the initial values.')
        print('There will be ', len(initial_values), ' simulations run. Will run parallel in sets of ',
              MAX_PROCESSES, '.', sep='')

        init = datetime.datetime.now()

        with Pool(MAX_PROCESSES) as p:
            tot_steps = sum(p.map(runAll, initial_values))

        end_t = datetime.datetime.now()

        if len(initial_values) > MAX_FILES:
            zipFiles(flsInDir(SUB_DIR), delete_after=False)
            deleteFolder()

        print('Did ', len(initial_values), ' simulations, with a total of  ', tot_steps,
              ' steps in ', end_t - init, '.', sep='')
    except FileNotFoundError:
        runAll()

    print('\nDone!')
