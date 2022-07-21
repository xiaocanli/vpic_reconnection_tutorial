"""Generate a XDMF meta file for VPIC fields and hydro data
"""
import numpy as np
import os
import sys
from lxml import etree

# data_format = "Binary"
data_format = "HDF"
compound_data = True
endian = "Little"  # Depending on machine

def get_vpic_info():
    """Get information of the VPIC simulation
    """
    with open('./info') as f:
        content = f.readlines()
    f.close()
    vpic_info = {}
    for line in content[1:]:
        if "=" in line:
            line_splits = line.split("=")
        elif ":" in line:
            line_splits = line.split(":")

        tail = line_splits[1].split("\n")
        vpic_info[line_splits[0].strip()] = float(tail[0])
    return vpic_info


def get_tframe_max(fields_interval):
    """Get maximum time frame for field and hydro dump
    """
    fdirs = [x[0] for x in os.walk("field_hdf5")]
    tindices = [int(fdir[13:]) for fdir in fdirs[1:]]
    tindices = np.asarray(tindices)
    return int(tindices.max() // fields_interval)


def scalar_field(field, species, dims, tindex, grid2):
    """Write for scalar field
    """
    var_name = field["name"]
    att2 = etree.SubElement(grid2, "Attribute",
                            attrib={'Name': var_name,
                                    'AttributeType': "Scalar",
                                    'Center': "Node"})
    dataitem = etree.SubElement(att2, "DataItem",
                                attrib={'Format': data_format,
                                        'ItemType': "Uniform",
                                        'DataType': "Float",
                                        'Precision': "4",
                                        'Endian': endian,
                                        'Dimensions': dims})
    var_component = field["vars"][0]
    if data_format == "Binary":
        dataitem.text = "data/" + var_component + "_" + str(tindex) + ".gda"
    else:
        if var_component in ["cbx", "cby", "cbz", "ex", "ey", "ez"]:
            fdir = "field_hdf5/"
            fname = "fields_" + str(tindex) + ".h5"
        else:
            fdir = "hydro_hdf5/"
            if species in ["e", "electron"]:
                sname = "electron"
            else:
                sname = "ion"
            fname = "hydro_" + sname + "_" + str(tindex) + ".h5"
        dataitem.text = (fdir + "T." + str(tindex) + "/" + fname +
                         ":/Timestep_" + str(tindex) + "/" + var_component)


def vector_field(field, species, dims, tindex, grid2):
    """Write for vector field
    """
    var_name = field["name"]
    dims_vector = dims + " 3"
    att2 = etree.SubElement(grid2, "Attribute",
                            attrib={'Name': var_name,
                                    'AttributeType': "Vector",
                                    'Center': "Node"})
    vec = etree.SubElement(att2, "DataItem",
                           attrib={'ItemType': "Function",
                                   'Dimensions': dims_vector,
                                   'Function': "JOIN($0, $1, $2)"})
    if var_name in ["E", "B"]:
        fdir = "field_hdf5/"
        fname = "fields_" + str(tindex) + ".h5"
    else:
        fdir = "hydro_hdf5/"
        if species in ["e", "electron"]:
            sname = "electron"
        else:
            sname = "ion"
        fname = "hydro_" + sname + "_" + str(tindex) + ".h5"
    for var_component in field["vars"]:
        dataitem = etree.SubElement(vec, "DataItem",
                                    attrib={'Format': data_format,
                                            'ItemType': "Uniform",
                                            'DataType': "Float",
                                            'Precision': "4",
                                            'Endian': endian,
                                            'Dimensions': dims})
        if data_format == 'Binary':
            dataitem.text = "data/" + var_component + "_" + str(tindex) + ".gda"
        else:
            dataitem.text = (fdir + "T." + str(tindex) + "/" + fname +
                             ":/Timestep_" + str(tindex) + "/" + var_component)


def tensor_field(field, species, dims, tindex, grid2):
    """Write for tensor field
    """
    var_name = field["name"]
    dims_tensor6 = dims + " 6"
    dims_tensor9 = dims + " 9"
    if len(field["vars"]) == 6:
        atype = "Tensor6"
        dims_data = dims_tensor6
        fun = "JOIN($0 $1 $2 $3 $4 $5)"
    elif len(field["vars"]) == 9:
        atype = "Tensor"
        dims_data = dims_tensor9
        fun = "JOIN($0 $1 $2 $3 $4 $5 $6 $7 $8)"
    att2 = etree.SubElement(grid2, "Attribute",
                            attrib={'Name': var_name,
                                    'AttributeType': atype,
                                    'Center': "Node"})
    vec = etree.SubElement(att2, "DataItem",
                           attrib={'ItemType': "Function",
                                   'Dimensions': dims_data,
                                   'Function': fun})
    fdir = "hydro_hdf5/"
    if species in ["e", "electron"]:
        sname = "electron"
    else:
        sname = "ion"
    fname = "hydro_" + sname + "_" + str(tindex) + ".h5"
    for var_component in field["vars"]:
        dataitem = etree.SubElement(vec, "DataItem",
                                    attrib={'Format': data_format,
                                            'ItemType': "Uniform",
                                            'DataType': "Float",
                                            'Precision': "4",
                                            'Endian': endian,
                                            'Dimensions': dims})
        if data_format == 'Binary':
            dataitem.text = "data/" + var_component + "_" + str(tindex) + ".gda"
        else:
            dataitem.text = (fdir + "T." + str(tindex) + "/" + fname +
                             ":/Timestep_" + str(tindex) + "/" + var_component)


def fields_list_binary():
    """Name list for binary format data
    """
    if compound_data:
        fields_list = [{"name": "E", "vars": ["ex", "ey", "ez"]},
                       {"name": "B", "vars": ["bx", "by", "bz"]},
                       {"name": "j", "vars": ["jx", "jy", "jz"]},
                       {"name": "ve", "vars": ["vex", "vey", "vez"]},
                       {"name": "ue", "vars": ["uex", "uey", "uez"]},
                       {"name": "vi", "vars": ["vix", "viy", "viz"]},
                       {"name": "ui", "vars": ["uix", "uiy", "uiz"]},
                       {"name": "pe", "vars": ["pe-xx", "pe-xy", "pe-xz",
                                               "pe-yx", "pe-yy", "pe-yz",
                                               "pe-zx", "pe-zy", "pe-zz"]},
                       {"name": "pi", "vars": ["pi-xx", "pi-xy", "pi-xz",
                                               "pi-yx", "pi-yy", "pi-yz",
                                               "pi-zx", "pi-zy", "pi-zz"]},
                       {"name": "asbJ", "vars": ["absJ"]},
                       {"name": "ne", "vars": ["ne"]},
                       {"name": "ni", "vars": ["ni"]}]
    else:
        fields_list = [{"name": "bx", "vars": ["bx"]},
                       {"name": "by", "vars": ["by"]},
                       {"name": "bz", "vars": ["bz"]},
                       {"name": "ex", "vars": ["ex"]},
                       {"name": "ey", "vars": ["ey"]},
                       {"name": "ez", "vars": ["ez"]},
                       {"name": "jx", "vars": ["jx"]},
                       {"name": "jy", "vars": ["jy"]},
                       {"name": "jz", "vars": ["jz"]},
                       {"name": "uex", "vars": ["uex"]},
                       {"name": "uey", "vars": ["uey"]},
                       {"name": "uez", "vars": ["uez"]},
                       {"name": "uix", "vars": ["uix"]},
                       {"name": "uiy", "vars": ["uiy"]},
                       {"name": "uiz", "vars": ["uiz"]},
                       {"name": "vex", "vars": ["vex"]},
                       {"name": "vey", "vars": ["vey"]},
                       {"name": "vez", "vars": ["vez"]},
                       {"name": "vix", "vars": ["vix"]},
                       {"name": "viy", "vars": ["viy"]},
                       {"name": "viz", "vars": ["viz"]},
                       {"name": "absJ", "vars": ["absJ"]},
                       {"name": "ne", "vars": ["ne"]},
                       {"name": "ni", "vars": ["ni"]},
                       {"name": "vkappa", "vars": ["vkappa"]},
                       {"name": "ne_vkappa", "vars": ["ne_vkappa"]}]

    return fields_list


def fields_list_hdf5():
    """Name list for HDF5 format data
    """
    if compound_data:
        fields_list = [{"name": "B", "vars": ["cbx", "cby", "cbz"]},
                       {"name": "E", "vars": ["ex", "ey", "ez"]},
                       {"name": "je", "vars": ["jx", "jy", "jz"], "species": "e"},
                       {"name": "ji", "vars": ["jx", "jy", "jz"], "species": "i"},
                       {"name": "pe", "vars": ["px", "py", "pz"], "species": "e"},
                       {"name": "pi", "vars": ["px", "py", "pz"], "species": "i"},
                       {"name": "te", "vars": ["txx", "txy", "txz",
                                               "tyy", "tyz", "tzz"],
                        "species": "e"},
                       {"name": "ti", "vars": ["txx", "txy", "txz",
                                               "tyy", "tyz", "tzz"],
                        "species": "i"},
                       {"name": "ke", "vars": ["ke"], "species": "e"},
                       {"name": "ki", "vars": ["ke"], "species": "i"},
                       {"name": "ne", "vars": ["rho"], "species": "e"},
                       {"name": "ni", "vars": ["rho"], "species": "i"}]
    else:
        fields_list = [{"name": "cbx", "vars": ["cbx"]},
                       {"name": "cby", "vars": ["cby"]},
                       {"name": "cbz", "vars": ["cbz"]},
                       {"name": "ex", "vars": ["ex"]},
                       {"name": "ey", "vars": ["ey"]},
                       {"name": "ez", "vars": ["ez"]},
                       {"name": "jex", "vars": ["jx"], "species": "e"},
                       {"name": "jey", "vars": ["jy"], "species": "e"},
                       {"name": "jez", "vars": ["jz"], "species": "e"},
                       {"name": "jix", "vars": ["jx"], "species": "i"},
                       {"name": "jiy", "vars": ["jy"], "species": "i"},
                       {"name": "jiz", "vars": ["jz"], "species": "i"},
                       {"name": "pex", "vars": ["px"], "species": "e"},
                       {"name": "pey", "vars": ["py"], "species": "e"},
                       {"name": "pez", "vars": ["pz"], "species": "e"},
                       {"name": "pix", "vars": ["px"], "species": "i"},
                       {"name": "piy", "vars": ["py"], "species": "i"},
                       {"name": "piz", "vars": ["pz"], "species": "i"},
                       {"name": "ke", "vars": ["ke"], "species": 'e'},
                       {"name": "ki", "vars": ["ke"], "species": 'i'},
                       {"name": "ne", "vars": ["rho"], "species": 'e'},
                       {"name": "ni", "vars": ["rho"], "species": 'i'}]

    return fields_list


def main():
    """business logic for when running this module as the primary one!"""
    vpic_info = get_vpic_info()
    if data_format == "Binary":
        dims = (str(int(vpic_info["nz"])) + " " +
                str(int(vpic_info["ny"])) + " " +
                str(int(vpic_info["nx"])))
        origin = (str(-0.5*vpic_info["Lz/de"]) + " " +
                  str(-0.5*vpic_info["Ly/de"]) + " 0")
        dxdydz = (str(vpic_info["dz/de"]) + " " +
                  str(vpic_info["dy/de"]) + " " +
                  str(vpic_info["dx/de"]))
    else:
        dims = (str(int(vpic_info["nx"])) + " " +
                str(int(vpic_info["ny"])) + " " +
                str(int(vpic_info["nz"])))
        origin = ("0 " +
                  str(-0.5*vpic_info["Ly/de"]) + " " +
                  str(-0.5*vpic_info["Lz/de"]))
        dxdydz = (str(vpic_info["dx/de"]) + " " +
                  str(vpic_info["dy/de"]) + " " +
                  str(vpic_info["dz/de"]))
    fields_interval = vpic_info["fields_interval"]
    tstart = 0
    tend = get_tframe_max(fields_interval)
    nframes = tend - tstart + 1
    tfields = vpic_info["dt*wpe"] * fields_interval # in 1/wpe
    tpoints = ""
    for i in range(tstart, tend+1):
        tpoints += str(tfields * i) + " "

    if data_format == 'Binary':
        fields_list = fields_list_binary()
    else:
        fields_list = fields_list_hdf5()

    xmlns_xi = {'xi': "http://www.w3.org/2001/XInclude"}
    root = etree.Element("Xdmf", nsmap=xmlns_xi, Version="2.0")

    domain = etree.SubElement(root, "Domain")

    topology = etree.SubElement(domain, "Topology",
                                attrib={'name': "topo",
                                        'TopologyType': "3DCoRectMesh",
                                        'Dimensions': dims})
    geometry = etree.SubElement(domain, "Geometry",
                                attrib={'name': "geo",
                                        'Type': "ORIGIN_DXDYDZ"})
    geometry.append(etree.Comment(" Origin "))
    dataitem1 = etree.SubElement(geometry, "DataItem",
                                 attrib={'Format': "XML",
                                         'Dimensions': "3"})
    dataitem1.text = origin
    geometry.append(etree.Comment(" DxDyDz "))
    dataitem2 = etree.SubElement(geometry, "DataItem",
                                 attrib={'Format': "XML",
                                         'Dimensions': "3"})
    dataitem2.text = dxdydz
    grid = etree.SubElement(domain, "Grid",
                            attrib={'Name': "TimeSeries",
                                    'GridType': "Collection",
                                    'CollectionType': "Temporal"})
    time = etree.SubElement(grid, "Time", TimeType="HyperSlab")
    dataitem = etree.SubElement(time, "DataItem",
                                attrib={'Format': "XML",
                                        'NumberType': "Float",
                                        'Dimensions': str(nframes)})
    dataitem.text = tpoints

    for tframe in range(tstart, tend+1):
        grid2 = etree.SubElement(grid, "Grid",
                                 attrib={'Name': "T" + str(tframe),
                                         'GridType': "Uniform"})
        topo2 = etree.SubElement(grid2, "Topology", Reference="/Xdmf/Domain/Topology[1]")
        geo2 = etree.SubElement(grid2, "Geometry", Reference="/Xdmf/Domain/Geometry[1]")
        tindex = int(tframe * fields_interval)
        for field in fields_list:
            if "species" in field:
                species = field["species"]
            else:
                species = "field"
            if len(field["vars"]) == 1:
                scalar_field(field, species, dims, tindex, grid2)
            if len(field["vars"]) == 3:
                vector_field(field, species, dims, tindex, grid2)
            if len(field["vars"]) > 3:
                tensor_field(field, species, dims, tindex, grid2)

    header = '<?xml version="1.0"?>\n'
    header += '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'

    vpic_xmf = etree.tostring(root, pretty_print=True, encoding='unicode')
    vpic_xmf = header + vpic_xmf
    with open('./vpic.xdmf', 'wb') as f:
        f.write(vpic_xmf.encode("utf-8"))

if __name__ == "__main__":
    main()
