/*****************************************************************************\
* This source file is part of mpisee
* Copyright (C) 2024  Ioannis Vardas - vardas@par.tuwien.ac.at
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>
******************************************************************************/
#include <cstdint>
#include <iostream>
#include <sqlite3.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <vector>
#include "utils.h"
#include "create_db.h"
#include <fcntl.h>


int countTable(sqlite3 *db, const std::string& tablename) {
    sqlite3_stmt* stmt;
    int count = 0;
    std::string checkSql = "SELECT COUNT(*) FROM " + tablename;
    sqlite3_prepare_v2(db, checkSql.c_str(), -1, &stmt, nullptr);

    if (sqlite3_step(stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(stmt, 0);
    }
    sqlite3_finalize(stmt);
    return count;
}


void printDataDetails(sqlite3* db) {
    sqlite3_stmt* stmt;

    std::string sql = "SELECT d.rank, m.machine, c.name, c.size, o.operation, "
                      "d.buffer_size_min, d.buffer_size_max, d.calls, d.time "
                      "FROM data d "
                      "JOIN mappings m ON d.rank = m.id "
                      "JOIN comms c ON d.comm_id = c.id "
                      "JOIN operations o ON d.operation_id = o.id";

    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int rank = sqlite3_column_int(stmt, 0);
            const char* machine = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
            const char* commName = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
            const char* commSize = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3));
            const char* operation = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 4));
            int bufferSizeMin = sqlite3_column_int(stmt, 5);
            int bufferSizeMax = sqlite3_column_int(stmt, 6);
            int calls = sqlite3_column_int(stmt, 7);
            double time = sqlite3_column_double(stmt, 8);

            std::cout << "Rank: " << rank
                      << ", Machine: " << (machine ? machine : "NULL")
                      << ", Comm Name: " << (commName ? commName : "NULL")
                      << ", Comm Size: " << (commSize ? commSize : "NULL")
                      << ", MPI Operation: " << (operation ? operation : "NULL")
                      << ", Buffer Size: " << bufferSizeMin << " - " << bufferSizeMax
                      << ", Calls: " << calls
                      << ", Time: " << time << std::endl;
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
}


void printData(sqlite3* db) {
    sqlite3_stmt* stmt;

    std::string sql = "SELECT c.name, c.size, d.rank, o.operation, "
                      "d.buffer_size_min, d.buffer_size_max, d.calls, d.time "
                      "FROM data d "
                      "JOIN comms c ON d.comm_id = c.id "
                      "JOIN operations o ON d.operation_id = o.id "
                      "ORDER BY c.name";  // Order by Comm Name

    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        // Print header
        std::cout << std::left << std::setw(15) << "Comm Name" << std::setw(15) << "Comm Size"
                  << std::setw(10) << "Rank" << std::setw(20) << "Operation"
                  << std::setw(25) << "Buffer Size Range" << std::setw(15) << "Calls"
                  << std::setw(20) << "Time" << std::endl;

        // Print rows
        while (sqlite3_step(stmt) == SQLITE_ROW) {
          std::stringstream bufferSizeStream;
            bufferSizeStream << sqlite3_column_int(stmt, 4) << " - " << sqlite3_column_int(stmt, 5);
            std::string bufferSize = bufferSizeStream.str();

            std::cout << std::left << std::setw(15) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0))
                      << std::setw(15) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1))
                      << std::setw(10) << sqlite3_column_int(stmt, 2)
                      << std::setw(20) << reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3))
                      << std::setw(25) << bufferSize //sqlite3_column_int(stmt, 4) << " - " << sqlite3_column_int(stmt, 5)
                      << std::setw(15) << sqlite3_column_int(stmt, 6)
                      << std::setw(20) << sqlite3_column_double(stmt, 7) << std::endl;
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
}



// Function to execute SQL command
void executeSQL(sqlite3* db, const std::string& sql, const char* name) {
    char* errMsg = nullptr;
    int rc = sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &errMsg);

    if (rc != SQLITE_OK) {
        std::cerr << "SQL error: " << errMsg << std::endl;
        sqlite3_free(errMsg);
    }
    // } else {
    //     std::cout << name << " successfully" << std::endl;
    // }
}


int getCommId(sqlite3* db, const std::string& commName) {
    sqlite3_stmt* stmt;
    int commId = -1;  // Default to an invalid ID
    // Prepare SQL query to select the ID from comms where name and size match
    std::string sql = "SELECT id FROM comms WHERE name = '" + commName + "'";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_text(stmt, 1, commName.c_str(), -1, SQLITE_STATIC);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            commId = sqlite3_column_int(stmt, 0);  // Get the id from the query result
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
    return commId; // return the retrieved ID
}

int getMappingId(sqlite3* db, const std::string& machineName) {
    sqlite3_stmt* stmt;
    int mappingId = -1;  // Default to an invalid ID

    std::string sql = "SELECT id FROM mappings WHERE machine = '" + machineName + "'";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_text(stmt, 1, machineName.c_str(), -1, SQLITE_STATIC);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            mappingId = sqlite3_column_int(stmt, 0);  // Get the id from the query result
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
    return mappingId;
}

void setMetadata(sqlite3 *db, const std::string &key,
                 const std::string &value) {


  std::string sql = "INSERT INTO metadata (key, value) VALUES ('" + key +
                    "', '" + value + "')";

    executeSQL(db, sql, "Metadata");
}

void insertMetadata(sqlite3 *db, char *mpi_lib, int size,
                    char *cmd[MAX_ARGS], int ac, int mpisee_major_v, int mpisee_minor_v,
                    char *build_date, char *build_time, const char *env) {
  std::string lib = mpi_lib;
  std::string sz = std::to_string(size);
  std::string mpisee_version =
      std::to_string(mpisee_major_v) + "." + std::to_string(mpisee_minor_v);

  std::string mpisee_date =
      std::string(build_date) + ", " + std::string(build_time);

  if( env != NULL )
    std::string env_var = env;
  std::string combinedString;
  for (int i = 0; i < ac && i< MAX_ARGS && cmd[i] != NULL; ++i) {
    if (i > 0) {
      combinedString += " ";
    }
    combinedString += cmd[i];
  }

  auto now = std::chrono::system_clock::now();

  // Convert time_point to time_t for easier manipulation
  std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);

  // Convert time_t to string representation
  std::string profile_date = std::ctime(&now_time_t);

  setMetadata(db, "MPI Library", lib);
  setMetadata(db, "Processes", sz);
  setMetadata(db, "Run command", combinedString);
  setMetadata(db, "mpisee version", mpisee_version);
  setMetadata(db, "mpisee build date", mpisee_date);
  setMetadata(db, "Profile date", profile_date);

}


void printMetadata(sqlite3* db) {
    sqlite3_stmt* stmt;
    std::string sql = "SELECT key, value FROM metadata";
    std::cout << "Printing Metadata\n";
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        const char *key =
            reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0));

        const char *value =
            reinterpret_cast<const char *>(sqlite3_column_text(stmt, 1));
        std::cout << key << ": " << value << '\n';
      }
    }
    sqlite3_finalize(stmt);
}

void createTables(sqlite3* db) {

    // Create Mapping Table
    const char* MappingTable =
        "CREATE TABLE IF NOT EXISTS mappings ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "machine TEXT);";
    executeSQL(db, MappingTable, "Mappings Table created");

    // Create Execution Time Table
    const char* ExecTimeTable =
        "CREATE TABLE IF NOT EXISTS exectimes ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "time REAL);";
    executeSQL(db, ExecTimeTable, "Execution Time Table created");

    // Create Metadata Table
    const char* Metadata =
        "CREATE TABLE IF NOT EXISTS metadata ("
        "key TEXT PRIMARY KEY, "
        "value TEXT);";
    executeSQL(db, Metadata, "Metadata Table created");

    // Create MPI Operations Table
    const char* MPIOpsTable =
        "CREATE TABLE IF NOT EXISTS operations ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "operation TEXT);";
    executeSQL(db, MPIOpsTable, "MPI Operations Table created ");

    // Create Comms Table
    const char* CommsTable =
        "CREATE TABLE IF NOT EXISTS comms ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "name TEXT UNIQUE, "
        "size INTEGER);";
    executeSQL(db, CommsTable, "Communicator Table created");

    // Create MPI Operations Volume Table
    // const char* MPIOpsVolumeTable =
    //     "CREATE TABLE IF NOT EXISTS ops_volume ("
    //     "id INTEGER PRIMARY KEY AUTOINCREMENT, "
    //     "operation_id INTEGER, "
    //     "rank INTEGER,"
    //     "comm_id INTEGER, "
    //     "volume INTEGER,"
    //     "FOREIGN KEY (operation_id) REFERENCES operations (id), "
    //     "FOREIGN KEY (comm_id) REFERENCES comms (id),"
    //     "FOREIGN KEY (rank) REFERENCES mappings (id));";
    // executeSQL(db, MPIOpsVolumeTable, "MPI Operations Volume Table created");

    // Create Data Table
    const char* DataTable =
        "CREATE TABLE IF NOT EXISTS data ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "rank INTEGER, "
        "comm_id INTEGER, "
        "operation_id INTEGER, "
        "buffer_size_min INTEGER, "
        "buffer_size_max INTEGER, "
        "calls INTEGER, "
        "time REAL, "
        "volume INTEGER, "
        "FOREIGN KEY (operation_id) REFERENCES operations (id), "
        "FOREIGN KEY (comm_id) REFERENCES comms (id), "
        "FOREIGN KEY (rank) REFERENCES mappings (id));";
    executeSQL(db, DataTable, "Data Table created");
}

// Function to insert into mappings
void insertIntoMappings(sqlite3 *db, const std::string &machine) {
  std::string insertSql;
  insertSql = "INSERT INTO mappings (id, machine) VALUES (0, '" + machine + "')";
  executeSQL(db, insertSql, "INSERT INTO mappings");
}

void BatchInsertIntoMappings(sqlite3 *db, const std::vector<std::string>& machines) {
  // Start a transaction
  executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

  for (const auto& machine : machines) {
    std::string insertSql;
    insertSql = "INSERT INTO mappings (machine) VALUES ('" + machine + "')";
    executeSQL(db, insertSql, "INSERT INTO mappings");
  }

  // Commit the transaction
  executeSQL(db, "END TRANSACTION", "End Transaction");
}

// Function to insert into mappings
void insertIntoTimes(sqlite3 *db, const double time) {
  std::string insertSql;
  insertSql = "INSERT INTO exectimes (id, time) VALUES (0, '" + std::to_string(time) + "')";
  executeSQL(db, insertSql, "INSERT INTO exectimes");
}

void BatchInsertIntoTimes(sqlite3 *db, const std::vector<double> times) {
  // Start a transaction
  executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

  for (const auto& time : times) {
    std::string insertSql;
    insertSql = "INSERT INTO exectimes (time) VALUES ('" + std::to_string(time) + "')";
    executeSQL(db, insertSql, "INSERT INTO exectimes");
  }

  // Commit the transaction
  executeSQL(db, "END TRANSACTION", "End Transaction");
}

int insertIntoComms(sqlite3 *db, const std::string &name, int size ) {
    // Insert or ignore based on unique name
    int commId = 0;
    sqlite3_stmt *insertStmt;
    sqlite3_stmt *getIdStmt;
    std::string insertSql = "INSERT OR IGNORE INTO comms (name, size) VALUES (?, ?)";
    std::string getIdSql = "SELECT id FROM comms WHERE name = ?";
    sqlite3_prepare_v2(db, insertSql.c_str(), -1, &insertStmt, nullptr);
    sqlite3_bind_text(insertStmt, 1, name.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int(insertStmt, 2, size);
    sqlite3_step(insertStmt);
    sqlite3_finalize(insertStmt);

    // Get the ID of the comm
    sqlite3_prepare_v2(db, getIdSql.c_str(), -1, &getIdStmt, nullptr);
    sqlite3_bind_text(getIdStmt, 1, name.c_str(), -1, SQLITE_STATIC);
    if (sqlite3_step(getIdStmt) == SQLITE_ROW) {
        commId = sqlite3_column_int(getIdStmt, 0);
    }
    sqlite3_finalize(getIdStmt);
    return commId;
}

std::vector<int> CommsInsert(sqlite3 *db, const std::vector<CommData>& comms) {
   sqlite3_stmt *insertStmt, *getIdStmt;
   std::vector<int> ids;
   int rc;
   std::string insertSql =
     "INSERT OR IGNORE INTO comms (name, size) VALUES (?, ?)";
   std::string getIdSql = "SELECT id FROM comms WHERE name = ?";

   rc = sqlite3_prepare_v2(db, insertSql.c_str(), -1, &insertStmt, nullptr);
   if (rc != SQLITE_OK) {
     std::cerr << "Failed to prepare insert statement: " << sqlite3_errmsg(db)
               << std::endl;

     return ids; // Return an empty vector or handle the error as appropriate
   }
   sqlite3_prepare_v2(db, getIdSql.c_str(), -1, &getIdStmt, nullptr);
   if (rc != SQLITE_OK) {
     std::cerr << "Failed to prepare getId statement: " << sqlite3_errmsg(db)
               << std::endl;
             sqlite3_finalize(insertStmt); // Clean up
     return ids; // Return an empty vector or handle the error as appropriate
   }

   // Start transaction
   executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

   for (const auto& comm : comms) {
     sqlite3_bind_text(insertStmt, 1, comm.name.c_str(), -1, SQLITE_STATIC);
     sqlite3_bind_int(insertStmt, 2, comm.size);
     rc = sqlite3_step(insertStmt);
     if (rc != SQLITE_DONE) {
       std::cerr << "Insert operation failed: " << sqlite3_errmsg(db)
                 << std::endl;
       break;
     }
     sqlite3_reset(insertStmt); // Reset the statement to insert next record


     // Get the ID of the communicator
     sqlite3_bind_text(getIdStmt, 1, comm.name.c_str(), -1, SQLITE_STATIC);
     if (sqlite3_step(getIdStmt) == SQLITE_ROW) {
       int id = sqlite3_column_int(getIdStmt, 0);
       ids.push_back(id); // Store the ID
     } else {
       std::cerr << "Failed to get ID: " << sqlite3_errmsg(db) << std::endl;
     }
     sqlite3_reset(getIdStmt);
   }

   // Finalize statement and commit transaction
   sqlite3_finalize(insertStmt);
   sqlite3_finalize(getIdStmt);
   executeSQL(db, "END TRANSACTION", "End Transaction");
   return ids;
}



// Functions to insert into operations
void insertIntoOperations(sqlite3* db, const std::string& operation) {
  int count;
  count = countTable(db, "operations");
  std::string insertSql;
  if (count == 0) {
    insertSql = "INSERT INTO operations (id, operation) VALUES (0, '" + operation + "')";
  } else {
    insertSql = "INSERT INTO operations (operation) VALUES ('" + operation + "')";

  }
  executeSQL(db, insertSql, "INSERT INTO operations");
}

void insertIntoOperationsEmpty(sqlite3 *db, const std::string &operation) {
  std::string insertSql;
  insertSql = "INSERT INTO operations (id, operation) VALUES (0, '" + operation + "')";
  executeSQL(db, insertSql, "INSERT INTO operations");
}

void BatchInsertIntoOperations(sqlite3 *db,
                               const std::vector<std::string> &operations) {
  // Start a transaction
  executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

  for (const auto& operation : operations) {
    std::string insertSql;
      // The table is not empty, let SQLite auto-increment the id
    insertSql = "INSERT INTO operations (operation) VALUES ('" + operation + "')";
    executeSQL(db, insertSql, "INSERT INTO operations");
  }

  // Commit the transaction
  executeSQL(db, "END TRANSACTION", "End Transaction");
}

// Functions to insert into data
void insertIntoData(sqlite3* db, int rank, int commId, int operationId,
                    int bufferSizeMax, int bufferSizeMin, int calls, double time) {
  std::string insertSql;
  insertSql = "INSERT INTO data (rank, comm_id, operation_id, buffer_size_max, buffer_size_min, calls, time) VALUES ("
                      + std::to_string(rank) + ", " + std::to_string(commId) + ", "
                      + std::to_string(operationId) + ", " + std::to_string(bufferSizeMax) + ", "
                      + std::to_string(bufferSizeMin) + ", " + std::to_string(calls) + ", " + std::to_string(time) + ")";
  executeSQL(db, insertSql, "INSERT INTO data");
}

void insertIntoDataEntry(std::vector<DataEntry> &entries, int rank, int commId,
                         int operationId, int bufferSizeMin, int bufferSizeMax,
                         int calls, double time, uint64_t volume) {
  // Create a new DataEntry object and add it to the vector
  // Use fields to initialize the object
    DataEntry entry = {
    .rank = rank,
    .commId = commId,
    .operationId = operationId,
    .bufferSizeMin = bufferSizeMin,
    .bufferSizeMax = bufferSizeMax,
    .calls = calls,
    .time = time,
    .volume = volume
  };
  // DataEntry entry = {rank, commId, operationId, bufferSizeMin, bufferSizeMax, calls, time};
    entries.push_back(entry);
}

void insertIntoVolEntry(std::vector<VolEntry> &entries, int operationId,
                        int rank, int commId, uint64_t volume) {
    VolEntry entry = {
    .operationId = operationId,
    .rank = rank,
    .commId = commId,
    .volume = volume
  };
    entries.push_back(entry);
}

void executeBatchInsert(sqlite3* db, const std::vector<DataEntry>& entries) {

    const std::string insertSql = "INSERT INTO data (rank, comm_id, operation_id, buffer_size_min, buffer_size_max, calls, time, volume) VALUES (?, ?, ?, ?, ?, ?, ?, ?)";

    sqlite3_stmt *stmt;
    int result = sqlite3_prepare_v2(db, insertSql.c_str(), -1, &stmt, nullptr);
    if (result != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        return;
    }

    // Start a single transaction
    executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

    for (const auto& entry : entries) {
        sqlite3_bind_int(stmt, 1, entry.rank);
        sqlite3_bind_int(stmt, 2, entry.commId);
        sqlite3_bind_int(stmt, 3, entry.operationId);
        sqlite3_bind_int(stmt, 4, entry.bufferSizeMin);
        sqlite3_bind_int(stmt, 5, entry.bufferSizeMax);
        sqlite3_bind_int(stmt, 6, entry.calls);
        sqlite3_bind_double(stmt, 7, entry.time);
        sqlite3_bind_int64(stmt, 8, entry.volume);

        result = sqlite3_step(stmt);
        if (result != SQLITE_DONE) {
            std::cerr << "Error during insert: " << sqlite3_errmsg(db) << std::endl;
        }

        sqlite3_reset(stmt); // Reset bindings for the next insertion
    }

    // Commit the transaction
    executeSQL(db, "END TRANSACTION", "End Transaction");

    sqlite3_finalize(stmt);

}


void printCommsTable(sqlite3* db) {
    sqlite3_stmt* stmt;

    // SQL query to select all records from the comms table
    std::string sql = "SELECT id, name, size FROM comms";

    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) == SQLITE_OK) {
        std::cout << "Contents of comms table:" << std::endl;
        std::cout << "ID\tName\tSize" << std::endl;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int id = sqlite3_column_int(stmt, 0);
            std::string name = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
            std::string size = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));

            std::cout << id << "\t" << name << "\t" << size << std::endl;
        }
        sqlite3_finalize(stmt);
    } else {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
}

void insertToVolume(sqlite3 *db, int operationId, int rank, int commId, uint64_t volume) {
    std::string insertSql;
    insertSql = "INSERT INTO ops_volume (operation_id, rank, comm_id, volume) VALUES ("
                + std::to_string(operationId) + ", " + std::to_string(rank) + ", "
                + std::to_string(commId) + ", " + std::to_string(volume) + ")";
    executeSQL(db, insertSql, "INSERT INTO ops_volume");
}

void batchInsertToVolume(sqlite3 *db, const std::vector<VolEntry> &entries) {
    const std::string insertSql = "INSERT INTO ops_volume (operation_id, rank, comm_id, volume) VALUES (?, ?, ?, ?)";

    sqlite3_stmt *stmt;
    int result = sqlite3_prepare_v2(db, insertSql.c_str(), -1, &stmt, nullptr);
    if (result != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        return;
    }

    // Start a single transaction
    executeSQL(db, "BEGIN TRANSACTION", "Start Transaction");

    for (const auto& entry : entries) {
        sqlite3_bind_int(stmt, 1, entry.operationId);
        sqlite3_bind_int(stmt, 2, entry.rank);
        sqlite3_bind_int(stmt, 3, entry.commId);
        sqlite3_bind_int(stmt, 4, entry.volume);

        result = sqlite3_step(stmt);
        if (result != SQLITE_DONE) {
            std::cerr << "Error during insert: " << sqlite3_errmsg(db) << std::endl;
        }

        sqlite3_reset(stmt); // Reset bindings for the next insertion
    }

    // Commit the transaction
    executeSQL(db, "END TRANSACTION", "End Transaction");

    sqlite3_finalize(stmt);
}



// Function to get the current time and date, create a stringstream to
// hold the output, create the condensed string and return the filename
std::string createFilename(std::string prefix, std::string suffix)
{
    // Get the current time
    std::time_t time = std::time(nullptr);
    std::tm tm = *std::localtime(&time);

    // Create a stringstream to hold the output
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d%H%M%S");

    // Prepend _mpisee and append .db to the string
    std::string filename = prefix + oss.str() + suffix;

    // Get the condensed string and return it
    return filename;
}


// Function to open SQLite database exclusively with retry mechanism
sqlite3* openSQLiteDBExclusively(std::string prefix, std::string suffix,
                                 int maxRetries, std::string &dbname)
{
    int retryCount = 0;
    std::string filename;
    sqlite3* db = NULL;
    int fd;

    do {
        filename = createFilename(prefix,suffix);
        // Use this file lock to ensure that the file is not opened by another process
        fd = open(filename.c_str(), O_RDWR | O_CREAT | O_EXCL, 0664);
        if (fd != -1) {
            int rc = sqlite3_open_v2(filename.c_str(), &db,
                                     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE |
                                     SQLITE_OPEN_FULLMUTEX, NULL);
            if (rc == SQLITE_OK) {
                std::cout << "mpisee: Opened database: " << filename <<
                  " exclusively after " << retryCount << " retries" << std::endl;
                dbname = filename;
                close(fd); // The handle (db) is now owned by SQLite
                return db;
            } else {
                mcpt_abort("Error opening database: %s",sqlite3_errmsg(db));
                if (db) {
                    close(fd);
                    sqlite3_close(db);
                    db = NULL;
                }
            }
            break;
        } else if (errno != EEXIST) {
            mcpt_abort("Error creating file: %s", filename.c_str());
            return NULL;
        }
        // File exists, sleep for a second, and retry
        sleep(1);
        ++retryCount;
    } while (retryCount < maxRetries);

    return NULL;
}
