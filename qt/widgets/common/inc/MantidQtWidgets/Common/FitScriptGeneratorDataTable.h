// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidQtWidgets/Common/IndexTypes.h"

#include <string>
#include <vector>

#include <QEvent>
#include <QModelIndex>
#include <QObject>
#include <QPainter>
#include <QPersistentModelIndex>
#include <QString>
#include <QStyleOptionViewItem>
#include <QStyledItemDelegate>
#include <QTableWidget>
#include <QTableWidgetItem>

namespace MantidQt {
namespace MantidWidgets {

/**
 * This class represents the table widget which holds domain data for the
 * FitScriptGenerator interface. This table has four columns:
 * Workspace Name, Workspace Index, Start X, End X.
 *
 * This table has been manually created and derived from QTableWidget to allow
 * the table rows to be highlighted when a hover event occurs.
 */
class FitScriptGeneratorDataTable : public QTableWidget {
  Q_OBJECT

public:
  enum ColumnIndex {
    WorkspaceName = 0,
    WorkspaceIndex = 1,
    StartX = 2,
    EndX = 3
  };

  FitScriptGeneratorDataTable(QWidget *parent = nullptr);
  ~FitScriptGeneratorDataTable() = default;

  std::string workspaceName(FitDomainIndex row) const;
  MantidWidgets::WorkspaceIndex workspaceIndex(FitDomainIndex row) const;
  double startX(FitDomainIndex row) const;
  double endX(FitDomainIndex row) const;

  std::vector<FitDomainIndex> selectedRows() const;

  void removeDomain(std::string const &workspaceName,
                    MantidWidgets::WorkspaceIndex workspaceIndex);
  void addDomain(QString const &workspaceName,
                 MantidWidgets::WorkspaceIndex workspaceIndex, double startX,
                 double endX);

  void formatSelection();
  void resetSelection();

signals:
  void itemExited(int newRowIndex);

private slots:
  void handleItemClicked(QTableWidgetItem *item);

private:
  bool eventFilter(QObject *widget, QEvent *event) override;
  QPersistentModelIndex hoveredRowIndex(QEvent *event);

  int indexOfDomain(std::string const &workspaceName,
                    MantidWidgets::WorkspaceIndex workspaceIndex) const;

  QString getText(FitDomainIndex row, int column) const;

  void setSelectedXValue(double xValue);

  int m_selectedRow;
  int m_selectedColumn;
  double m_selectedValue;
  QPersistentModelIndex m_lastIndex;
};

/**
 * This class is used for formating the type of data allowed in each of the
 * tables columns. It is also used for setting various column properties, and
 * will paint a row when it is hovered over.
 */
class CustomItemDelegate : public QStyledItemDelegate {
  Q_OBJECT

public:
  enum class DelegateType { Double, Int, String };

  CustomItemDelegate(FitScriptGeneratorDataTable *parent = nullptr,
                     DelegateType const &type = DelegateType::Double);

private slots:
  void handleItemEntered(QTableWidgetItem *item);
  void handleItemExited(int newRowIndex);

private:
  QWidget *createEditor(QWidget *parent, QStyleOptionViewItem const &option,
                        QModelIndex const &index) const override;
  void paint(QPainter *painter, QStyleOptionViewItem const &option,
             QModelIndex const &index) const override;

  FitScriptGeneratorDataTable *m_tableWidget;
  int m_hoveredIndex;
  DelegateType m_type;
};

} // namespace MantidWidgets
} // namespace MantidQt
