import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { MotifMatchStatusComponent } from './motif-match-status.component';

describe('MotifMatchStatusComponent', () => {
  let component: MotifMatchStatusComponent;
  let fixture: ComponentFixture<MotifMatchStatusComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ MotifMatchStatusComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(MotifMatchStatusComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
