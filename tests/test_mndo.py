
from fitter import mndo

def test_params():

    parameters = {}
    parameters["O"] = {}
    parameters["O"]["USS"] = 666.0
    filename = "_tmp_test_params"

    mndo.set_params(parameters)

    atoms = [
    'O',
    'N',
    'C',
    'N',
    'N',
    'H',]

    coords = [
        [ -0.0593325887, 1.2684201211  , 0.0095178503  ],
        [ 1.1946293712 , 1.771776509   , 0.0001229152  ],
        [ 1.9590217387 , 0.7210517427  , -0.0128069641 ],
        [ 1.2270421979 , -0.4479406483 , -0.0121559722 ],
        [ 0.0119302176 , -0.1246338293 , 0.0012973737  ],
        [ 3.0355546734 , 0.7552313348  , -0.0229864829 ],
    ]

    inptxt = mndo.get_input(atoms, coords, 0, "testfile")

    f = open(filename, 'w')
    f.write(inptxt)
    f.write(inptxt)
    f.close()

    calculations = mndo.run_mndo_file(filename)

    for lines in calculations:
        properties = mndo.get_properties(lines)

        idx = get_index(lines, "USS")
        line = stdout[idx]
        line = line.split()

        value = float(line[-1])

        assert value == 666.0

    return

