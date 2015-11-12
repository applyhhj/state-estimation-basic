package thu.instcloud.app.se.mpdata;

import thu.instcloud.app.se.common.Constants;

import java.util.ArrayList;
import java.util.List;

import static thu.instcloud.app.se.common.Utils.Common.readStringFromFile;

/**
 * Created on 2015/11/6.
 */
public class MPData {

//    private static Logger logger = LoggerFactory.getLogger(MPData.class);

    private BranchData branchData;

    private BusData busData;

    private GeneratorData genData;

    private double sbase;

    private String version;

    private int nBranch;

    private int nBus;

    public MPData(String filepath) {

        branchData = new BranchData();

        busData = new BusData();

        genData = new GeneratorData();

        importData(filepath);

        busData.reorderBusNumbers(branchData);

        nBranch = branchData.getN();

        nBus = busData.getN();

    }

    private void importData(String filepath) {

        List<String> fileContent = readStringFromFile(filepath);

        String entry;

        List<String> data = new ArrayList<String>();

        boolean end;

        int i = 0;

        while (i < fileContent.size()) {

            entry = fileContent.get(i);

            if (entry.startsWith("#")) {

                i++;

                continue;

            }

            end = false;

            if (entry.contains(Constants.MPC.MPC_VERSION)) {

                version = fileContent.get(++i);

                if (fileContent.get(++i).contains(Constants.MPC.MPC_VERSION_END)) {

                    end = true;

                }

            } else if (entry.contains(Constants.MPC.MPC_BASEMVA)) {

                sbase = Double.parseDouble(fileContent.get(++i));

                if (fileContent.get(++i).contains(Constants.MPC.MPC_BASEMVA_END)) {

                    end = true;

                }

            } else if (entry.contains(Constants.MPC.MPC_BUS)) {

                data.clear();

                while (++i < fileContent.size()) {

                    if (fileContent.get(i).contains(Constants.MPC.MPC_BUS_END)) {

                        end = true;

                        break;

                    } else {

                        data.add(fileContent.get(i));

                    }

                }

                if (end) {

                    if (!busData.loadData(data)) {

                        break;

                    }

                }

            } else if (entry.contains(Constants.MPC.MPC_GEN)) {

                data.clear();

                while (++i < fileContent.size()) {

                    if (fileContent.get(i).contains(Constants.MPC.MPC_GEN_END)) {

                        end = true;

                        break;

                    } else {

                        data.add(fileContent.get(i));

                    }

                }

                if (end) {

                    if (!genData.loadData(data)) {

                        break;

                    }

                }

            } else if (entry.contains(Constants.MPC.MPC_BRANCH)) {

                data.clear();

                while (++i < fileContent.size()) {

                    if (fileContent.get(i).contains(Constants.MPC.MPC_BRANCH_END)) {

                        end = true;

                        break;

                    } else {

                        data.add(fileContent.get(i));

                    }

                }

                if (end) {

                    if (!branchData.loadData(data)) {

                        break;

                    }

                }

            } else {

                end = true;

            }

            if (!end) {

                System.out.print("Error: Section " + entry + " has no end");

                return;

            }

            i++;

        }

    }

    public double getSbase() {
        return sbase;
    }

    public BranchData getBranchData() {
        return branchData;
    }

    public BusData getBusData() {
        return busData;
    }

    public String getVersion() {
        return version;
    }

    public int getnBranch() {
        return nBranch;
    }

    public int getnBus() {
        return nBus;
    }
}

