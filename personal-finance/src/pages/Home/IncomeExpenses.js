import React, { Component } from 'react';

import Assets from '../../components/Assets';
import Income from '../../components/Income';
import Expenses from '../../components/Expenses';
import Storage from '../../components/Storage';

import CONFIG from '../../config';

class IncomeExpenses extends Component {
  constructor(props) {
    super(props);

    this.state = {
      all_items: Storage.getItems(null, props.items),
      assets: Storage.getItems(CONFIG.ITEM_TYPE.ASSET, props.items),
      liabilities: Storage.getItems(CONFIG.ITEM_TYPE.LIABILITY, props.items),
      income: Storage.getItems(CONFIG.ITEM_TYPE.INCOME, props.items),
      expenses: Storage.getItems(CONFIG.ITEM_TYPE.EXPENSE, props.items),
    };
  }

  render() {
    return (
      <React.Fragment>
        <div className="container">
          <br />

          <div className="row">
            <div className="col-sm">
              <Income items={this.state.income} />
            </div>

            <div className="col-sm">
              <Expenses items={this.state.expenses} />
            </div>
          </div>
        </div>
      </React.Fragment>
    );
  }
}

export default IncomeExpenses;
